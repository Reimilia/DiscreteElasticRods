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
	restPos.resize(3, nNodes);
	for (int i = 0; i < nNodes; i++)
	{
		nodes.push_back(particles[i]);
		restPos.col(i) = particles[i].pos;
	}

	rods.clear();
	rods.shrink_to_fit();
	restLength.resize(nNodes - 1);
	for (int i = 0; i < nNodes - 1; i++)
	{
		double l = (nodes[i + 1].pos - nodes[i].pos).norm();
		rods.push_back(ElasticRodSegment(i,i+1,0,l));
		restLength[i] = l;
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
		stencils[i].restlength = l;
	}

	params = para;

	beta = 1.0;

	//Precompute something

	// The first one material frame need to set up manually
	// The best choice is the body template coordinate frame at the first centerline.
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

	// Compute the twist via quasi static assumption
	updateMaterialCurvature();
	restCurvature.resize(3, 2*(nNodes - 2));
	for (int i = 0; i < nNodes - 2; i++)
	{
		restCurvature.col(2 * i) = stencils[i].prevCurvature;
		restCurvature.col(2 * i + 1) = stencils[i].nextCurvature;
	}
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

bool ElasticRod::buildConfiguration(Eigen::VectorXd & pos, Eigen::VectorXd & vel)
{
	int nparticles = (int)nodes.size();
	pos.resize(3 * nparticles);
	vel.resize(3 * nparticles);
	for (int i = 0; i < nparticles; i++)
	{
		pos.segment<3>(3 * i) = nodes[i].pos;
		vel.segment<3>(3 * i) = nodes[i].vel;
	}
	return true;
}

bool ElasticRod::unbuildConfiguration(const Eigen::VectorXd & pos, const Eigen::VectorXd & vel)
{
	int nparticles = (int)nodes.size();
	for (int i = 0; i < nparticles; i++)
	{
		nodes[i].prevpos = pos;
		nodes[i].pos = pos.segment<3>(3 * i);
		nodes[i].vel = vel.segment<3>(3 * i);
	}
	return true;
}

double ElasticRod::computeTotalEnergy()
{
	// Twist + Bend
	// Bend is what we realized in World of Goo Mile Stone
	double energy = 0.0;
	int nstencils = (int)stencils.size();
	for (int i = 0; i < nstencils; i++)
	{
		double e1 = (stencils[i].prevCurvature - restCurvature.col(2 * i)).dot(rods[i].bendingModulus * (stencils[i].prevCurvature - restCurvature.col(2 * i)));
		double e2 = (stencils[i].nextCurvature - restCurvature.col(2 * i + 1)).dot(rods[i + 1].bendingModulus * (stencils[i].prevCurvature - restCurvature.col(2 * i + 1)));
		energy += (e1 + e2) / 2 / stencils[i].restlength;
		energy += beta * (rods[i + 1].theta - rods[i].theta) * (rods[i + 1].theta - rods[i].theta) / stencils[i].restlength;
	}
	return energy;
}

void ElasticRod::computeCenterlineForces(Eigen::VectorXd &f)
{
	/*
	Compute -dE as the force

	The difficulty here is boundary conditions
	*/
	Eigen::Matrix2d J;
	J << 0, -1, 1, 0;

	int nstencils = (int)stencils.size();
	int nparticles = (int)nodes.size();
	double dEdtheta = (stencils[nstencils - 1].nextCurvature.dot(J * rods[nstencils].bendingModulus* (stencils[nstencils - 1].nextCurvature - restCurvature.col(2 * nstencils - 1))) +
		2 * beta * (rods[nstencils].theta - rods[nstencils - 1].theta)) / stencils[nstencils - 1].restlength;

	f.resize(3 * nparticles);

	for (int k = 0; k < nstencils; k++)
	{
		Eigen::Vector2d coef1, coef2;
		coef1 = rods[k].bendingModulus * (stencils[k].prevCurvature - restCurvature.col(2*k))/ stencils[k].restlength;
		coef2 = rods[k + 1].bendingModulus * (stencils[k].nextCurvature - restCurvature.col(2 * k + 1))/ stencils[k].restlength;

		Eigen::Vector3d e1, e2;
		Eigen::Matrix3d dkb1, dkb2;
		e1 = nodes[k + 1].pos - nodes[k].pos;
		e2 = nodes[k + 2].pos - nodes[k + 1].pos;
		dkb1 = (2 * VectorMath::crossProductMatrix(e2) + stencils[k].kb * e2.transpose()) / (restLength[k]*restLength[k+1] + e1.dot(e2));
		dkb2 = (2 * VectorMath::crossProductMatrix(e1) + stencils[k].kb * e1.transpose()) / (restLength[k] * restLength[k + 1] + e1.dot(e2));
		
		Eigen::MatrixXd psi;
		psi.resize(nparticles, 3);

		for (int i = k; i < nparticles; i++)
		{
			//compute domega
			
			// j=k
			// j=k+1
		}
	}

	for (int i = 0; i < nparticles; i++)
	{
		Eigen::Vector3d psi(0,0,0);
		if (i > 1) psi -= stencils[i - 2].kb / 2 / rods[i - 1].length;
		if (i > 0 && i <= nstencils) psi += (stencils[i-1].kb / 2 / rods[i].length - stencils[i-1].kb/2/rods[i-1].length);
		if (i < nstencils) psi += stencils[i].kb / 2 / rods[i].length;
		f.segment<3>(3 * i) += dEdtheta * psi;
	}
}

void ElasticRod::computeEnergyThetaDifferentialAndHessian(Eigen::VectorXd & dE, Eigen::VectorXd & lowerH, Eigen::VectorXd & centerH, Eigen::VectorXd & upperH)
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
	dE.setZero();
	lowerH.setZero();
	centerH.setZero();
	upperH.setZero();
	
	Eigen::Matrix2d J;
	J << 0, -1, 1, 0;

	for (int i = 0; i < nrods; i++)
	{
		dE[i] = lowerH[i] = centerH[i] = upperH[i] = 0;
		if (i < nrods-1)
		{
			dE[i] -= 2 * beta * (rods[i + 1].theta - rods[i].theta) / stencils[i].restlength;
			double dw = (stencils[i].prevCurvature).dot(J*rods[i].bendingModulus* (stencils[i].prevCurvature - restCurvature.col(2 * i)));
			dE[i] +=  dw / stencils[i].restlength;
			upperH[i] = -2 * beta / stencils[i].restlength;
			centerH[i] += 2 * beta / stencils[i].restlength + (stencils[i].prevCurvature.dot(J.transpose()*rods[i].bendingModulus*J*stencils[i].prevCurvature) + dw) / stencils[i].restlength;
		}
		if (i > 0)
		{
			dE[i] += 2 * beta * (rods[i].theta - rods[i - 1].theta) / stencils[i - 1].restlength;
			double dw = (stencils[i - 1].nextCurvature).dot(J*rods[i].bendingModulus* (stencils[i - 1].nextCurvature - restCurvature.col(2 * i - 1)));
			dE[i] +=  dw / stencils[i - 1].restlength;
			lowerH[i] = -2 * beta / stencils[i - 1].restlength;
			centerH[i] += 2 * beta / stencils[i - 1].restlength + (stencils[i - 1].nextCurvature.dot(J.transpose()*rods[i].bendingModulus*J*stencils[i - 1].nextCurvature) + dw) / stencils[i - 1].restlength;
		}
		
	}
	
}


void ElasticRod::updateAfterTimeIntegration()
{
	/*
	Updates the following:
	1. rods' length (Although we have inextensible constraint)
	2. bishop frame
	3. check quasi static frame constraint and update material frame(i.e. theta)
	
	*/

	int nrods = (int)rods.size();
	for (int i = 0; i < nrods; i++)
	{
		rods[i].length = (nodes[i + 1].pos - nodes[i].pos).norm();
		if (i > 0)
		{
			stencils[i - 1].length = rods[i].length + rods[i - 1].length;
		}
	}
	updateBishopFrame();
	updateQuasiStaticFrame();
	updateMaterialCurvature();
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
		computeEnergyThetaDifferentialAndHessian(dE, lowerH, centerH, upperH);
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
