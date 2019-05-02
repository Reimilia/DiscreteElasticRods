#include "ElasticRod.h"


ElasticRod::ElasticRod()
{
}

ElasticRod::ElasticRod(std::vector<ElasticRod> rods, std::vector<Particle, Eigen::aligned_allocator<Particle>> nodes)
{
	//Create bending stencils 
	
	//Precompute something
}

ElasticRod::~ElasticRod()
{
}

const Eigen::Vector3d ElasticRod::renderPos(Eigen::Vector3d & point, int index)
{
	/*
	For rendering purpose:
	Need to know the actual position in world coordinate. 
	The rigid body template transformation need to be combined with material curvature in centerline.
	*/
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

	// This is bad, but not that important
	int maxNewtonIter = 10;
	double newtonTolerance = 1e-6;

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


	for (int t = 0; t < maxNewtonIter; t++)
	{
		computeEnergyDifferentialAndHessian(dE, lowerH, centerH, upperH);
		if (dE.norm() < newtonTolerance)
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
	// The first one need to set up manually
	t0 = nodes[1].pos - nodes[0].pos;
	t0 = t0 / t0.norm();
	// u0 $\perp$ t0
	u0 = Eigen::Vector3d(t0[2]-t0[1],t0[0]-t0[2],t0[1]-t0[0]);
	u0 = u0 / u0.norm();
	v0 = u0.cross(t0);

	rods[0].t = t0;
	rods[0].u= u0;
	rods[0].v = v0;


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
