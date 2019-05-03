#include "ElasticRod.h"


ElasticRod::ElasticRod()
{
	nodes.clear();
	nodes.shrink_to_fit();
	rods.clear();
	rods.shrink_to_fit();
	stencils.clear();
	stencils.shrink_to_fit();
}

ElasticRod::ElasticRod(std::vector<Particle, Eigen::aligned_allocator<Particle>> &particles, SimParameters para, BoundaryCondition bc, std::vector<double> &initialTheta)
{
	//Create rods and bending stencils 
	int nNodes = (int)particles.size();
	nodes.clear();
	nodes.shrink_to_fit();
	for (int i = 0; i < nNodes; i++)
	{
		nodes.push_back(particles[i]);
	}

	rods.clear();
	rods.shrink_to_fit();
	for (int i = 0; i < nNodes - 1; i++)
	{
		double l = (nodes[i + 1].pos - nodes[i].pos).norm();
		rods.push_back(ElasticRodSegment(i,i+1,0,l));
		rods[i].theta = initialTheta[i];
	}

	stencils.clear();
	stencils.shrink_to_fit();
	for (int i = 0; i < nNodes - 2; i++)
	{
		Eigen::Vector3d e1, e2;
		e1 = nodes[i + 1].pos - nodes[i].pos;
		e2 = nodes[i + 2].pos - nodes[i + 1].pos;
		Eigen::Vector3d kb = 2 * e1.cross(e2) / (e1.norm()*e2.norm() + e1.dot(e2));
		double l = rods[i].length + rods[i + 1].length;
		stencils.push_back(BendingStencil(i, i + 1, i + 2, kb));
		stencils[i].length = l;
	}

	params = para;

	//Precompute something

	// Store rest position and velocity

	// The first one material frame need to set up manually
	t0 = nodes[1].pos - nodes[0].pos;
	t0 = t0 / t0.norm();
	// u0 $\perp$ t0
	u0 = Eigen::Vector3d(t0[2] - t0[1], t0[0] - t0[2], t0[1] - t0[0]);
	u0 = u0 / u0.norm();
	v0 = u0.cross(t0);

	rods[0].t = t0;
	rods[0].u = u0;
	rods[0].v = v0;

	bcStats = bc;

	// Compute the rest material frame
	updateBishopFrame();
	
	// Compute Material Curvature
	updateQuasiStaticFrame();

	// Compute the rest material curvature
	updateMaterialCurvature();
}

ElasticRod::~ElasticRod()
{
}

const Eigen::Vector3d ElasticRod::renderPos(Eigen::Vector3d & point, int index)
{
	/*
	For rendering purpose:
	Need to know the actual position in world coordinate. 
	The rigid body template transformation need to be combined with material frame in centerline.
	Basically, one lazy way is to compute interpolated theta based on the position, 
	and then rotate along t-axis via the interpolated theta value

	The input of this function is the rigid body point position before apply twist!
	We do not bind any structure of centerline in the class.
	*/
	
	Eigen::Vector3d pos;
	int p1;
	p1 = rods[index].p1;
	double dist = (point - nodes[p1].pos).dot(rods[index].t);

	if (dist<0 || dist>rods[index].length)
		return point;

	return VectorMath::rotationMatrix(rods[index].t * rods[index].theta / rods[index].length * dist) * point;
}

void ElasticRod::computeEnergyDifferentialAndHessian(Eigen::VectorXd & dE, Eigen::VectorXd & lowerH, Eigen::VectorXd & centerH, Eigen::VectorXd & upperH)
{
	/*
	Only for inner d.o.f. (theta)
	All d.o.f. that need to be computed via time integrator will be put in outer scope.
	*/
	int nrods = (int)rods.size();

	dE.resize(nrods);
	lowerH.resize(nrods);
	centerH.resize(nrods);
	upperH.resize(nrods);
	
	Eigen::Matrix2d J;
	J << 0, -1, 1, 0;

	for (int i = 0; i < nrods; i++)
	{
		dE[i] = lowerH[i] = centerH[i] = upperH[i] = 0;
		if (i < nrods-1)
		{
			dE[i] -= 2 * beta * (rods[i + 1].theta - rods[i].theta) / stencils[i].length;
			double dw = stencils[i].prevCurvature.dot(J*rods[i].bendingModulus* (stencils[i].prevCurvature - restCurvature.row(2 * i)));
			dE[i] +=  dw / stencils[i].length;
			upperH[i] = -2 * beta / stencils[i].length;
			centerH[i] += 2 * beta / stencils[i].length + (stencils[i].prevCurvature.dot(J.transpose()*rods[i].bendingModulus*J* stencils[i].prevCurvature) + dw) / stencils[i].length;
		}
		if (i > 0)
		{
			dE[i] += 2 * beta * (rods[i].theta - rods[i - 1].theta) / stencils[i - 1].length;
			double dw = stencils[i - 1].nextCurvature.dot(J*rods[i].bendingModulus* (stencils[i - 1].nextCurvature - restCurvature.row(2 * i - 1)));
			dE[i] +=  dw / stencils[i - 1].length;
			lowerH[i] = -2 * beta / stencils[i - 1].length;
			centerH[i] += 2 * beta / stencils[i - 1].length + (stencils[i - 1].nextCurvature.dot(J.transpose()*rods[i].bendingModulus*J*stencils[i - 1].nextCurvature) + dw) / stencils[i - 1].length;
		}
		
	}
	
}

void ElasticRod::updateQuasiStaticFrame()
{
	/*
	Using Newton's method to solve non-linear equation to determine twist angle theta in each rod element
	Note that it shall incoporates boundary conditions
	*/

	// This is bad, but in case we need this
	//int NewtonMaxIter = 10;
	//double NewtonTolerance = 1e-6;

	int nrods = (int)rods.size();

	Eigen::VectorXd dE;
	Eigen::VectorXd upperH, lowerH, centerH;

	// Build configuration
	Eigen::VectorXd thetas;
	thetas.resize(nrods);
	for (int i = 0; i < nrods; i++)
	{
		thetas[i] = rods[i].theta;
	}


	for (int t = 0; t < params.NewtonMaxIters; t++)
	{
		computeEnergyDifferentialAndHessian(dE, lowerH, centerH, upperH);
		if (dE.norm() < params.NewtonTolerance)
		{
			break;
		}
		switch (bcStats)
		{
		case BC_FREE:
			break;
		case BC_FIXED_END:
			upperH[0] = 0;
			centerH[0] = centerH[nrods - 1] = 1.0;
			lowerH[nrods - 1] = 0;
		default:
			break;
		}
		MatrixMath::tridiagonalSolver(lowerH, centerH, upperH, dE);
		thetas -= dE;
	}

	// Unbuild configuration
	for (int i = 0; i < nrods; i++)
	{
		rods[i].theta = thetas[i];
	}
}

void ElasticRod::updateBishopFrame()
{
	/*	
	Compute Discrete Parallel Transportation to update Bishop Frame for centerline
	Move coordinate frame via the rotational matrix
	*/

	int nRods = (int)rods.size();


	// Now compute Bishop frame
	for (int i = 1; i < nRods; i++)
	{
		rods[i].t = nodes[i + 1].pos - nodes[i].pos;
		rods[i].t = (rods[i].t) / (rods[i].t).norm();
		Eigen::Vector3d n = (rods[i-1].t).cross(rods[i].t);

		rods[i].u = VectorMath::rotationMatrix(n*asin(n.norm()) / n.norm()) * rods[i - 1].u;
		rods[i].u = (rods[i].u) / (rods[i].u).norm();

		rods[i].v = (rods[i].t).cross(rods[i].u);
	}
}

void ElasticRod::updateMaterialCurvature()
{
	/*
	Compute material curvatures
	*/
	int nRods = (int)rods.size();
	for (int i = 0; i < nRods; ++i)
	{
		Eigen::Vector3d m1 = cos(rods[i].theta)* rods[i].u + sin(rods[i].theta)* rods[i].v;
		Eigen::Vector3d m2 = -sin(rods[i].theta)* rods[i].v + cos(rods[i].theta)* rods[i].u;
		if (i < nRods - 1)
		{
			stencils[i].prevCurvature = Eigen::Vector2d((stencils[i].kb).dot(m2), -(stencils[i].kb).dot(m1));
		}
		if (i > 0)
		{
			stencils[i - 1].nextCurvature = Eigen::Vector2d((stencils[i - 1].kb).dot(m2), -(stencils[i - 1].kb).dot(m1));
		}
		
	}
}
