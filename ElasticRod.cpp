#include "ElasticRod.h"
#include "SimParameters.h"
#include <iostream>

ElasticRod::ElasticRod()
{
	nodes.clear();
	nodes.shrink_to_fit();
	rods.clear();
	rods.shrink_to_fit();
	stencils.clear();
	stencils.shrink_to_fit();
}

ElasticRod::ElasticRod(std::vector<Particle, Eigen::aligned_allocator<Particle>> &particles, SimParameters para)
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
	if (nNodes == 1)
		return;

	rods.clear();
	rods.shrink_to_fit();
	restLength.resize(nNodes - 1);
	for (int i = 0; i < nNodes - 1; i++)
	{
		double l = (nodes[i + 1].pos - nodes[i].pos).norm();
		rods.push_back(ElasticRodSegment(i,i+1,params.rodDensity * l,l));
		rods[i].bendingModulus = params.rodBendingModulus;
		restLength[i] = l;
		restLength[i] = l;
		nodes[i].mass += params.rodDensity * l / 2;
		nodes[i + 1].mass += params.rodDensity * l / 2;
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
	
	leftRigidBody = rightRigidBody = -1;

	params = para;

	beta = params.rodTwistStiffness;

	//Precompute something

	// The first one material frame need to set up manually
	// The best choice is the body template coordinate frame at the first centerline.
	rods[0].t = nodes[1].pos - nodes[0].pos;
	rods[0].t = rods[0].t / rods[0].t.norm();
	rods[0].u = Eigen::Vector3d(rods[0].t[2] - rods[0].t[1], rods[0].t[0] - rods[0].t[2], rods[0].t[1] - rods[0].t[0]);
	rods[0].u = rods[0].u / rods[0].u.norm();
	rods[0].v = rods[0].t.cross(rods[0].u);
	rods[0].v /= rods[0].v.norm();

	switch (para.boundaryCondition)
	{
	case SimParameters::BC_FIXED_END:
		if (nNodes < 2)
			break;
		rods[0].theta = para.leftAngularVelocity * para.timeStep;

		rods[nNodes - 2].theta = para.rightAngularVelocity * para.timeStep;
	case SimParameters::BC_RIGIDBODY_END:
		// TODO: finish this
		break;
	case SimParameters::BC_FREE:
		break;
	default:
		break;
	}
	// Compute the bishop frame
	updateBishopFrame();

	// Compute the twist via quasi static assumption
	/*
	I do not fully understand why compute rest curvature goes first than quasi static frame
	*/
	updateMaterialCurvature();
	restCurvature.resize(2, 2*(nNodes - 2));
	for (int i = 0; i < nNodes - 2; i++)
	{
		restCurvature.col(2 * i) = stencils[i].prevCurvature;
		restCurvature.col(2 * i + 1) = stencils[i].nextCurvature;
	}

	// Compute Material Curvature
	updateQuasiStaticFrame();

	// update Material Frame
	updateMaterialCurvature();
}

ElasticRod::~ElasticRod()
{
}

bool ElasticRod::setSimulationParameters(SimParameters para)
{
	params = para;
	int nRods = (int)rods.size();
	// If clamped the bounary, need to recompute everything
	/*switch(para.boundaryCondition)
	{case SimParameters::BC_FIXED_END:
		if ((int)rods.size() < 1)
			break;
		rods[0].theta = 0;

		rods[nRods - 1].theta = 0;

		updateMaterialCurvature();
		for (int i = 0; i < nRods - 1; i++)
		{
			restCurvature.col(2 * i) = stencils[i].prevCurvature;
			restCurvature.col(2 * i + 1) = stencils[i].nextCurvature;
		}
		updateQuasiStaticFrame();

		updateMaterialCurvature();
		break;
	case SimParameters::BC_RIGIDBODY_END:
		// TODO: finish this
		break;
	case SimParameters::BC_FREE:
		break;
	default:
		break;
	}*/

	for (int i = 0; i < (int)rods.size(); i++)
	{
		rods[i].bendingModulus = params.rodBendingModulus;
	}
	
	return true;
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

	/*if (dist<0 || dist>rods[index].length)
	{
		return point;
	}*/

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
		nodes[i].prevpos = nodes[i].pos;
		nodes[i].pos = pos.segment<3>(3 * i);
		nodes[i].vel = vel.segment<3>(3 * i);
	}

	if (nparticles >= 3 && params.boundaryCondition == SimParameters::BC_FIXED_END)
	{
		rods[0].theta += params.timeStep * params.leftAngularVelocity;
		rods[nparticles-2].theta += params.timeStep * params.rightAngularVelocity;
	}
	
	updateAfterPosChange();
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
		double e2 = (stencils[i].nextCurvature - restCurvature.col(2 * i + 1)).dot(rods[i + 1].bendingModulus * (stencils[i].nextCurvature - restCurvature.col(2 * i + 1)));
		energy += (e1 + e2) / 2 / stencils[i].restlength;
		energy += beta * (rods[i + 1].theta - rods[i].theta) * (rods[i + 1].theta - rods[i].theta) / stencils[i].restlength;
	}
	std::cout << energy << std::endl;
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
	
	f.resize(3 * nparticles);
	f.setZero();

	/*Eigen::VectorXd dE, upperH, lowerH, centerH;
	Eigen::VectorXd thetas;
	thetas.resize(nstencils + 1);

	for (int i = 0; i < nstencils + 1; i++)
	{
		thetas[i] = rods[i].theta;
	}

	computeEnergyThetaDifferentialAndHessian(thetas, dE, lowerH, centerH, upperH);
	dE[0] = dE[nstencils] = 0;
	std::cout << "Sanity Check: dE/d\\theta=0 for all free rods : " << dE.norm() << std::endl;
	*/
	//std::cout << "1:" << f.norm() << std::endl;
	for (int k = 0; k < nstencils; k++)
	{
		Eigen::Vector2d coef1, coef2;
		coef1 = rods[k].bendingModulus * (stencils[k].prevCurvature - restCurvature.col(2 * k)) / stencils[k].restlength;
		coef2 = rods[k + 1].bendingModulus * (stencils[k].nextCurvature - restCurvature.col(2 * k + 1))/ stencils[k].restlength;


		Eigen::Vector3d e1, e2;
		Eigen::Matrix3d dkb1, dkb2;
		e1 = nodes[k + 1].pos - nodes[k].pos;
		e2 = nodes[k + 2].pos - nodes[k + 1].pos;
		dkb1 = (2 * VectorMath::crossProductMatrix(e2) + stencils[k].kb * e2.transpose()) / (restLength[k] * restLength[k + 1] + e1.dot(e2));
		dkb2 = (2 * VectorMath::crossProductMatrix(e1) - stencils[k].kb * e1.transpose()) / (restLength[k] * restLength[k + 1] + e1.dot(e2));
		
		//std::cout << dkb1 << std::endl;
	    //std::cout << dkb2 << std::endl;

		for (int i = 0; i <= k+2; i++)
		{
			//compute domega
			Eigen::Vector3d psi(0,0,0);
			Eigen::MatrixXd M;
			M.resize(2, 3);
			if (i >= nparticles) continue;
			// j=k
			M.row(1) = - cos(rods[k].theta)* rods[k].u - sin(rods[k].theta)* rods[k].v;
			M.row(0) = - sin(rods[k].theta)* rods[k].u + cos(rods[k].theta)* rods[k].v;
			// std::cout << M << std::endl;
			switch (i - k)
			{
			case 0:
				M = M * dkb1;
				break;
			case 1:
				M = M * (-dkb1 - dkb2);
				break;
			case 2:
				M = M * dkb2;
				break;
			default:
				M.setZero();
				break;
			}
			if (i >= 2 && i <= k + 1)
			{
				psi = psi - stencils[i - 2].kb / restLength[i - 1];
			}
			if (i >= 1 && i <= k)
			{
				psi = psi - stencils[i - 1].kb / restLength[i - 1] +  stencils[i - 1].kb / restLength[i];
			}

			if (i <= k-1)
			{
				psi = psi + stencils[i].kb / restLength[i];
			}
			psi /= 2;
			//Assemble force (watch out for minus sign!)
			f.segment<3>(3 * i) += (-(M - J * stencils[k].prevCurvature*psi.transpose()).transpose()*coef1);

			//j=k+1
			M.row(1) = - cos(rods[k+1].theta)* rods[k+1].u - sin(rods[k+1].theta)* rods[k+1].v;
			M.row(0) = - sin(rods[k+1].theta)* rods[k+1].u + cos(rods[k+1].theta)* rods[k+1].v;
			switch (i - k)
			{
			case 0:
				M = M * dkb1;
				break;
			case 1:
				M = M * (-dkb1 - dkb2);
				break;
			case 2:
				M = M * dkb2;
				break;
			default:
				M.setZero();
				break;
			}
			psi.setZero();
			if (i >= 2 && i <= k + 2)
			{
				psi = psi - stencils[i - 2].kb / restLength[i - 1];
			}
			if (i >= 1 && i <= k + 1) 
			{
				psi = psi - stencils[i - 1].kb / restLength[i - 1] + stencils[i - 1].kb / restLength[i];
			}
			if (i <= k)
			{
				psi = psi + stencils[i].kb / restLength[i];
			}
			psi /= 2;

			//Assemble force (watch out for minus sign!)
			f.segment<3>(3 * i) += (-(M - J * stencils[k].nextCurvature*psi.transpose()).transpose()*coef2);
			//std::cout << "2:" << f.norm() << std::endl;
		}


		
	}

	
	/* Test code, found numerical difference in force computation :(
	
	int i, j, k;
	i = 0; j = 0; k = 0;
	Eigen::MatrixXd W;
	
	auto func = [this](Eigen::MatrixXd &W,int i,int j,int k){
		W.resize(2, 3);
		W.setZero();

		Eigen::Matrix2d J;
		J << 0, -1, 1, 0;
		
		Eigen::Vector3d e1, e2;
		Eigen::Matrix3d dkb1, dkb2;
		e1 = nodes[k + 1].pos - nodes[k].pos;
		e2 = nodes[k + 2].pos - nodes[k + 1].pos;
		dkb1 = (2 * VectorMath::crossProductMatrix(e2) + stencils[k].kb * e2.transpose()) / (restLength[k] * restLength[k + 1] + e1.dot(e2));
		dkb2 = (2 * VectorMath::crossProductMatrix(e1) - stencils[k].kb * e1.transpose()) / (restLength[k] * restLength[k + 1] + e1.dot(e2));

		Eigen::Vector3d psi(0, 0, 0);
		Eigen::MatrixXd M;
		M.resize(2, 3);
		M.row(1) = - cos(rods[k].theta) * rods[k].u - sin(rods[k].theta) * rods[k].v;
		M.row(0) = - sin(rods[k].theta) * rods[k].u + cos(rods[k].theta) * rods[k].v;
		if (j == k)
		{
			switch (i - k)
			{
			case 0:
				M = M * dkb1;
				break;
			case 1:
				M = M * (-dkb1 - dkb2);
				break;
			case 2:
				M = M * dkb2;
				break;
			default:
				M.setZero();
				break;
			}
			if (i >= 2 && i <= k + 1)
			{
				psi = psi - stencils[i - 2].kb / restLength[i - 1];
			}
			std::cout << psi.transpose() << std::endl;
			if (i >= 1 && i <= k)
			{
				psi = psi - stencils[i - 1].kb / restLength[i - 1] + stencils[i - 1].kb / restLength[i];
			}
			std::cout << psi.transpose() << std::endl;
			if (i <= k - 1)
			{
				psi = psi + stencils[i].kb / restLength[i];
			}
			psi /= 2;
			std::cout << psi.transpose() << std::endl;
			std::cout << M << std::endl;
			W = (M - J * stencils[k].prevCurvature*psi.transpose());
		}
		else
		{
			M.row(1) = -cos(rods[k + 1].theta)* rods[k + 1].u - sin(rods[k + 1].theta)* rods[k + 1].v;
			M.row(0) = -sin(rods[k + 1].theta)* rods[k + 1].u + cos(rods[k + 1].theta)* rods[k + 1].v;
			switch (i - k)
			{
			case 0:
				M = M * dkb1;
				break;
			case 1:
				M = M * (-dkb1 - dkb2);
				break;
			case 2:
				M = M * dkb2;
				break;
			default:
				M.setZero();
				break;
			}
			psi.setZero();
			if (i >= 2 && i <= k + 2)
			{
				psi = psi - stencils[i - 2].kb / restLength[i - 1];
			}
			if (i >= 1 && i <= k + 1)
			{
				psi = psi - stencils[i - 1].kb / restLength[i - 1] + stencils[i - 1].kb / restLength[i];
			}
			if (i <= k)
			{
				psi = psi + stencils[i].kb / restLength[i];
			}
			psi /= 2;
			std::cout << psi.transpose() << std::endl;
			W = (M - J * stencils[k].nextCurvature*psi.transpose());
		}
		
	};

	Eigen::Vector2d w = stencils[0].prevCurvature;

	std::cout << "Original Curvature is : " << w.transpose() << std::endl;
	Eigen::VectorXd direction = Eigen::VectorXd::Random(3);
	direction.normalize();

	for (int k = 3; k <= 12; k++)
	{
		double eps = pow(10, -k);
		Eigen::MatrixXd dW;
		nodes[0].pos += direction*eps;
		int nrods = (int)rods.size();
		updateAfterPosChange();
		func(dW, 0, 0, 0);

		Eigen::Vector2d epsW = stencils[0].prevCurvature;
		std::cout << "Curvature is: " << epsW.transpose() << std::endl;
		std::cout << "EPS is: " << eps << std::endl;
		std::cout << "Norm of Finite Difference is: " << (epsW - w).norm() / eps << std::endl;
		std::cout << "Norm of Directinal Gradient is: " << (dW*direction).norm() << std::endl;
		std::cout << "The difference between above two is: " << ((epsW - w) / eps - dW * direction).norm() << std::endl << std::endl;

		nodes[0].pos -= direction * eps;
		updateAfterPosChange();
	
	}*/
	//std::cout << "3:" << f.norm() << std::endl;
	if (params.boundaryCondition == SimParameters::BC_FIXED_END)
	{
		double dEdtheta = 
			(stencils[nstencils - 1].nextCurvature.dot(J * rods[nstencils].bendingModulus* (stencils[nstencils - 1].nextCurvature - restCurvature.col(2 * nstencils - 1))) +
			2 * beta * (rods[nstencils].theta - rods[nstencils - 1].theta)) / stencils[nstencils - 1].restlength;

		std::cout << dEdtheta << std::endl;
		for (int i = 0; i < nparticles; i++)
		{
			Eigen::Vector3d psi;	
			psi.setZero();
			if (i > 1) psi = psi - stencils[i - 2].kb / restLength[i - 1];
			if (i > 0 && i <= nstencils) psi = psi + (stencils[i - 1].kb / restLength[i]) - (stencils[i - 1].kb / restLength[i - 1]);
			if (i < nstencils) psi = psi + (stencils[i].kb / restLength[i]);
			psi = psi / 2;
			f.segment<3>(3 * i) += (dEdtheta * psi);
		}

	}
	//std::cout << "4:" << f.norm() << std::endl;
}

void ElasticRod::computeEnergyThetaDifferentialAndHessian(const Eigen::VectorXd &theta, Eigen::VectorXd & dE, Eigen::VectorXd & lowerH, Eigen::VectorXd & centerH, Eigen::VectorXd & upperH)
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

	// Since theta changed, we need update material curvature!


	for (int i = 0; i < nrods; i++)
	{
		//std::cout << i << ": " << rods[i].theta << std::endl;
		dE[i] = lowerH[i] = centerH[i] = upperH[i] = 0;

		Eigen::Vector3d m1 = cos(theta[i])* rods[i].u + sin(theta[i])* rods[i].v;
		Eigen::Vector3d m2 = -sin(theta[i])* rods[i].u + cos(theta[i])* rods[i].v;
		if (i < nrods-1)
		{
			Eigen::Vector2d prevCurvature = Eigen::Vector2d((stencils[i].kb).dot(m2), -(stencils[i].kb).dot(m1));
			dE[i] -= 2 * beta * (theta[i + 1] - theta[i]) / stencils[i].length;
			double dw = prevCurvature.dot(J*rods[i].bendingModulus* (prevCurvature - restCurvature.col(2 * i)));
			dE[i] +=  dw / stencils[i].length;
			upperH[i] = -2 * beta / stencils[i].length;
			centerH[i] += 2 * beta / stencils[i].length +
				(prevCurvature.dot(J.transpose()*rods[i].bendingModulus*J*prevCurvature) - 
				(prevCurvature).dot(rods[i].bendingModulus* (prevCurvature - restCurvature.col(2 * i)))) / stencils[i].length;
		}
		if (i > 0)
		{
			Eigen::Vector2d nextCurvature = Eigen::Vector2d((stencils[i - 1].kb).dot(m2), -(stencils[i - 1].kb).dot(m1));
			dE[i] += 2 * beta * (theta[i] - theta[i - 1]) / stencils[i - 1].length;
			double dw = nextCurvature.dot(J*rods[i].bendingModulus* (nextCurvature - restCurvature.col(2 * i - 1)));
			dE[i] +=  dw / stencils[i - 1].length;
			lowerH[i] = -2 * beta / stencils[i - 1].length;
			centerH[i] += 2 * beta / stencils[i - 1].length +
				(nextCurvature.dot(J.transpose()*rods[i].bendingModulus*J*nextCurvature) - 
				nextCurvature.dot(rods[i].bendingModulus* (nextCurvature - restCurvature.col(2 * i - 1)))) / stencils[i - 1].length;
		}
		
	}


}


void ElasticRod::updateAfterPosChange()
{
	/*
	Updates the following:
	1. rods' length (Although we have inextensible constraint) and kb
	2. bishop frame
	3. check quasi static frame constraint and update material frame(i.e. theta)
	
	*/
	//std::cout << "Thetas:\n";
	int nrods = (int)rods.size();
	for (int i = 0; i < nrods; i++)
	{
		rods[i].length = (nodes[i + 1].pos - nodes[i].pos).norm();
		if (i > 0)
		{
			stencils[i - 1].length = rods[i].length + rods[i - 1].length;
		}
		//std::cout << rods[i].theta << std::endl;
	}
	for (int i = 0; i < nrods - 1; i++)
	{
		Eigen::Vector3d e1, e2;
		e1 = nodes[i + 1].pos - nodes[i].pos;
		e2 = nodes[i + 2].pos - nodes[i + 1].pos;
		stencils[i].kb = 2 * e1.cross(e2) / (restLength[i] * restLength[i+1] + e1.dot(e2));
		//std::cout << stencils[i].kb.transpose() << std::endl;
	}
	//compute t,u,v
	updateBishopFrame();
	//compute theta
	updateQuasiStaticFrame();
	//compute omega
	updateMaterialCurvature();
}

bool ElasticRod::updateQuasiStaticFrame()
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

	int t;
	for (t = 0; t < params.NewtonMaxIters; t++)
	{
		//std::cout << thetas.transpose() << std::endl;
		computeEnergyThetaDifferentialAndHessian(thetas, dE, lowerH, centerH, upperH);
		//std::cout << dE.norm() << std::endl;
		switch (params.boundaryCondition)
		{
		case SimParameters::BC_FREE:
			break;
		case SimParameters::BC_FIXED_END:
			/*
			Change matrix as:
			1 0 0 ...
			....
			... 0 0 1
			*/
			upperH[0] = 0;
			centerH[0] = centerH[nrods - 1] = 1.0;
			lowerH[nrods - 1] = 0;
			dE[0] = 0;
			dE[nrods - 1] = 0;
		default:
			break;
		}
		
		if (dE.norm() < params.NewtonTolerance)
		{
			std::cout << "Newton's Iteration converges in " << t << "steps.\n";
			//flag = true;
			break;
		}
		
		MatrixMath::tridiagonalSolver(lowerH, centerH, upperH, dE);
		thetas -= dE;
		

	}
	//std::cout << thetas.transpose() << std::endl;
	// Unbuild configuration
	if (t < params.NewtonMaxIters)
	{
		for (int i = 0; i < nrods; i++)
		{
			rods[i].theta = thetas[i];
		}
		updateMaterialCurvature();
		return true;
	}
	else
	{
		std::cout << "Update theta fails!\n";
		return false;
	}
}

void ElasticRod::updateBishopFrame()
{
	/*	
	Compute Discrete Parallel Transportation to update Bishop Frame for centerline
	Move coordinate frame via the rotational matrix
	*/

	int nRods = (int)rods.size();
	rods[0].t = nodes[1].pos - nodes[0].pos;
	rods[0].t = rods[0].t / rods[0].t.norm();
	rods[0].u = Eigen::Vector3d(rods[0].t[2] - rods[0].t[1], rods[0].t[0] - rods[0].t[2], rods[0].t[1] - rods[0].t[0]);
	rods[0].u = rods[0].u / rods[0].u.norm();
	rods[0].v = rods[0].t.cross(rods[0].u);
	rods[0].v /= rods[0].v.norm();
	

	// Now compute Bishop frame
	for (int i = 1; i < nRods; i++)
	{
		rods[i].t = nodes[i + 1].pos - nodes[i].pos;
		rods[i].t = (rods[i].t) / (rods[i].t).norm();
		Eigen::Vector3d n = (rods[i - 1].t).cross(rods[i].t);

		// Watchout!
		if (n.norm() < 1e-10)
		{
			rods[i].u = rods[i - 1].u;
		}
		else
		{
			if (rods[i].t.dot(rods[i - 1].t) > 0)
			{
				rods[i].u = VectorMath::rotationMatrix(n*asin(n.norm()) / n.norm()) * rods[i - 1].u;
			}
			else
			{
				rods[i].u = VectorMath::rotationMatrix(n*(M_PI - asin(n.norm())) / n.norm()) * rods[i - 1].u;
			}
			
		}
		//rods[i].u = VectorMath::rotationMatrix(stencils[i-1].kb) * rods[i-1].u;
		
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
		Eigen::Vector3d m2 = -sin(rods[i].theta)* rods[i].u + cos(rods[i].theta)* rods[i].v;
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
