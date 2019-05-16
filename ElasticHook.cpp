
#include "ElasticHook.h"
#include "TestModule.h"

#include <Eigen/SparseQR>



void ElasticHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen))
	{
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		/*if (ImGui::Button("Load", ImVec2((w - p) / 2.f, 0)))
		{
			if (!isPaused())
				pause();
			std::string filePath = igl::file_dialog_open();
			loadConfiguration(filePath);
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Save", ImVec2((w - p) / 2.f, 0)))
		{
			std::string filePath = igl::file_dialog_save();
			saveConfiguration(filePath);
		}*/

		//Test Part
		if (ImGui::Button("Test", ImVec2(-1, 0)))
		{
			testProcess();
		}
	}
	/*if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Combo("Click Adds", (int *)&params_.clickMode, "Particles\0\0");
		ImGui::Combo("Connector Type", (int *)&params_.connectorType, "Flexible Rods\0\0");
	}*/
	if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		//ImGui::Combo("Constraint Handling", (int *)&params_.constraintHandling, "Penalty Method\0Step and Project\0Lagrange Multipliers\0\0");
		//ImGui::Combo("Integrator", (int *)&params_.integrator, "Velocity Verlet\0Implicit Midpoint\0\0");
		ImGui::InputDouble("Timestep", &params_.timeStep);
		ImGui::InputDouble("Newton Tolerance", &params_.NewtonTolerance);
		ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
		//ImGui::InputDouble("Penalty Stiffness", &params_.penaltyStiffness);
	}
	if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
		ImGui::InputDouble("  Gravity g", &params_.gravityG);
		ImGui::Checkbox("Floor Enabled", &params_.floorEnabled);
		//ImGui::Checkbox("Bending Enabled", &params_.bendingEnabled);
		//ImGui::Checkbox("Twist Enabled", &params_.twistEnabled);
		//viewer.imgui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
	}


	if (ImGui::CollapsingHeader("New Particles", ImGuiTreeNodeFlags_DefaultOpen))
	{
		//ImGui::Checkbox("Is Fixed", &params_.particleFixed);
		ImGui::InputDouble("Mass", &params_.particleMass);
	}

	if (ImGui::CollapsingHeader("New Rods", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::InputDouble("Density", &params_.rodDensity);
		ImGui::InputDouble("Twisting Stiffness", &params_.rodTwistStiffness);
		ImGui::InputDouble("B(0,0)", &params_.rodBendingModulus(0, 0));
		ImGui::InputDouble("B(0,1)", &params_.rodBendingModulus(0, 1));
		ImGui::InputDouble("B(1,0)", &params_.rodBendingModulus(1, 0));
		ImGui::InputDouble("B(1,1)", &params_.rodBendingModulus(1, 1));
		ImGui::Combo("Boundary Conditions", (int *)&params_.boundaryCondition, "Free\0Clamped End\0\0");
	}

}


void ElasticHook::updateRenderGeometry()
{
	double baseradius = 0.01;
	double pulsefactor = 0.1;
	double pulsespeed = 50.0;

	//int sawteeth = 20;
	//double sawdepth = 0.1;
	//double sawangspeed = 10.0;

	double baselinewidth = 0.005;

	int numcirclewedges = 20;

	int totverts = 0;
	int totfaces = 0;

	// floor

	totverts += 5;
	totfaces += 4;

	for (RigidBodyInstance *rbi : bodies_)
	{
		totverts += rbi->getTemplate().getVerts().rows();
		totfaces += rbi->getTemplate().getFaces().rows();
	}

	for (ElasticRod *rod : rods_)
	{
		int nverts = rod->nodes.size();
		totverts += nverts * ballTemplate_->getVerts().rows() + (nverts - 1) * rodTemplate_->getVerts().rows();
		totfaces += nverts * ballTemplate_->getFaces().rows() + (nverts - 1) * rodTemplate_->getFaces().rows();
	}


	renderQ.resize(totverts, 3);
	renderF.resize(totfaces, 3);
	int voffset = 0;
	int foffset = 0;

	double eps = 1e-4;


	if (params_.floorEnabled)
	{
		/*for (int i = 0; i < 5; i++)
		{
		vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
		}*/

		renderQ.row(voffset++) =  Eigen::Vector3d(0, -1, 0);
		renderQ.row(voffset++) = Eigen::Vector3d(1e3, -1, 1e3);
		renderQ.row(voffset++) = Eigen::Vector3d(-1e3, -1, 1e3);
		renderQ.row(voffset++) = Eigen::Vector3d(-1e3, -1, -1e3);
		renderQ.row(voffset++) = Eigen::Vector3d(1e3, -1, -1e3);

		renderF.row(foffset) = Eigen::Vector3i(foffset + 0, foffset + 2, foffset + 1);
		renderF.row(foffset+1) = Eigen::Vector3i(foffset + 0, foffset + 3, foffset + 2);
		renderF.row(foffset+2) = Eigen::Vector3i(foffset + 0, foffset + 4, foffset + 3);
		renderF.row(foffset+3) = Eigen::Vector3i(foffset + 0, foffset + 1, foffset + 4);
		
		foffset += 4;
	}


	double basescale = 0.1;
	
	for (RigidBodyInstance *rbi : bodies_)
    {
        int nverts = rbi->getTemplate().getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (rbi->c + VectorMath::rotationMatrix(rbi->theta)*rbi->getTemplate().getVerts().row(i).transpose()).transpose();
        int nfaces = rbi->getTemplate().getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = rbi->getTemplate().getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }


	//Lazy rendering for rods
	//TODO: find a better way to render it
	Eigen::Vector3d dir0(1, 0, 0);
	for (ElasticRod *rod : rods_)
	{
		int nParticles = rod->nodes.size();
		for (int k = 0; k < nParticles; k++)
		{
			int nverts = ballTemplate_->getVerts().rows();
			for (int i = 0; i < nverts; i++)
				renderQ.row(voffset + i) = (rod->nodes[k].pos + ballTemplate_->getVerts().row(i).transpose()).transpose();
			int nfaces = ballTemplate_->getFaces().rows();
			for (int i = 0; i < nfaces; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					renderF(foffset + i, j) = ballTemplate_->getFaces()(i, j) + voffset;
				}
			}
			voffset += nverts;
			foffset += nfaces;
		}

		for (int k = 0; k < nParticles - 1; k++)
		{
			int nverts = rodTemplate_->getVerts().rows();
			ElasticRodSegment onebar = rod->getRodSegment(k);
			Eigen::Vector3d width(10 * onebar.length , 1, 1);
			Eigen::Vector3d c = (rod->getNodeSegment(k).pos + rod->getNodeSegment(k+1).pos) / 2;
			Eigen::Matrix3d rotation;
			rotation.col(0) = onebar.t;

			// TODO : Interpolate theta
			rotation.col(1) = cos(onebar.theta) * onebar.u + sin(onebar.theta) * onebar.v;
			rotation.col(2) = -sin(onebar.theta) * onebar.u + cos(onebar.theta) * onebar.v;
			//rotation.col(1) = onebar.u;
			//rotation.col(2) = onebar.v;
			Eigen::Vector3d theta = VectorMath::axisAngle(rotation);
			
			for (int i = 0; i < nverts; i++)
			{
				Eigen::Vector3d point = rodTemplate_->getVerts().row(i);
				double dist = (point - rod->getNodeSegment(k).pos).dot(onebar.t);
				Eigen::Vector3d pos;
				pos = c + VectorMath::rotationMatrix(theta)* point.cwiseProduct(width);
				//renderQ.row(voffset + i) = rod->renderPos(pos,k);
				renderQ.row(voffset + i) = pos;
			}	
			int nfaces = rodTemplate_->getFaces().rows();
			for (int i = 0; i < nfaces; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					renderF(foffset + i, j) = rodTemplate_->getFaces()(i, j) + voffset;
				}
			}
			voffset += nverts;
			foffset += nfaces;
		}

	}

}


void ElasticHook::initSimulation()
{
	time_ = 0;

	particles_.clear();
	particles_.shrink_to_fit();
	for (ElasticRod *it : rods_)
		delete it;
	rods_.clear();
	rods_.shrink_to_fit();
	isNewRod_ = false;
	loadMesh();
}

bool ElasticHook::mouseClicked(igl::opengl::glfw::Viewer & viewer, Eigen::Vector3d dir, int button)
{
	if (button != 2 || !isPaused())
		return false;

	message_mutex.lock();
	launch_ = true;
	Eigen::Matrix4f view = viewer.core.view;
	Eigen::Vector4f eye = view.inverse() * Eigen::Vector4f(0, 0, 0, 1.0f);
	for (int i = 0; i < 3; i++)
	{

		launchPos_[i] = eye[i] + dir[i];
		launchDir_[i] = dir[i];
	}

	message_mutex.unlock();
	return true;
}


void ElasticHook::tick()
{
	message_mutex.lock();
	{
		if (launch_ && isPaused())
		{
			addParticle(launchPos_[0], launchPos_[1], launchPos_[2]);
			launch_ = false;
		}
	}
	message_mutex.unlock();
	
}

bool ElasticHook::simulateOneStep()
{

	//Eigen::VectorXd q, v, prevq;
	//buildConfiguration(q, v);
	//bool isSucceed = numericalIntegration();
	//unbuildConfiguration(q, v);

	/*
	We need one time step of Lagrangian for constranints
	Then we need another time step for collision response.
	*/
	computeMassInverse(mInv_);
	bool isSucceed = numericalIntegration();
	time_ += params_.timeStep;
	return (false || !isSucceed);
}

void ElasticHook::addParticle(double x, double y, double z)
{
	Eigen::Vector3d newpos(x, y, z);
	Eigen::Vector3d zero(0, 0, 0);

	//RigidBodyInstance *rbi = new RigidBodyInstance(*ballTemplate_, newpos, zero, zero, zero, 1);
	//bodies_.push_back(rbi);
	if (isNewPause_)
	{
		std::cout << "Here!" << std::endl;
		particles_.clear();
		particles_.shrink_to_fit();
		isNewPause_ = false;
	}
	else
	{
		if ((int)rods_.size() > 0) 
		{
			delete rods_.back();
			rods_.pop_back();
		}
	}
	particles_.push_back(Particle(newpos, params_.particleMass,false, false));

	particles_[0].fixed = true;
	particles_[(int)particles_.size() - 1].fixed = true;
	rods_.push_back(new ElasticRod(particles_, params_));
	// Assign Boundary Condition
	// TODO: Assign Rigid body condition
	if (params_.boundaryCondition == SimParameters::BC_FIXED_END)
	{
		if ((int)particles_.size() > 1)
		{
			rods_.back()->assignClampedBoundaryCondition(0, M_PI);
		}
	}
	
	particles_[(int)particles_.size() - 1].fixed = false;
}

void ElasticHook::buildConfiguration(Eigen::VectorXd &pos, Eigen::VectorXd &vel)
{
	int nconfigurations = 0;
	for (int i = 0; i < (int)rods_.size(); i++)
	{
		nconfigurations += (3 * rods_[i]->nodes.size());
	}

	pos.resize(nconfigurations);
	vel.resize(nconfigurations);
	int offset = 0;
	for (int i = 0; i < (int)rods_.size(); i++)
	{
		Eigen::VectorXd q,v;
		rods_[i]->buildConfiguration(q,v);
		for (int k = 0; k < q.size(); k++)
		{
			pos(offset + k) = q(k);
			vel(offset + k) = v(k);
		}
		offset += q.size();
	}
}

void ElasticHook::unbuildConfiguration(const Eigen::VectorXd &pos, const Eigen::VectorXd &vel)
{
	int offset = 0;
	for (int i = 0; i < (int)rods_.size(); i++)
	{
		Eigen::VectorXd q, v;
		int ndofs = 3 * rods_[i]->nodes.size();
		q.resize(ndofs);
		v.resize(ndofs);
		for (int k = 0; k < ndofs; k++)
		{
			q(k) = pos(offset + k);
			v(k) = vel(offset + k);
		}
		rods_[i]->unbuildConfiguration(q, v);
		offset += ndofs;
	}
	
}

void ElasticHook::computeMassInverse(Eigen::SparseMatrix<double>& mInv)
{
	// First rods, then rigid bodies
	
	std::vector<Eigen::Triplet<double>> massInvTriplet;
	int offset = 0;

	for (int i = 0; i < (int)rods_.size(); i++)
	{
		ElasticRod & rod = *rods_[i];
		int nParticles = (int)rod.nodes.size();
		for (int k = 0; k < nParticles; k++)
		{
			double massInverse;
			if (rod.nodes[k].fixed)
				massInverse = 0;
			else
				massInverse = 1.0 / rod.nodes[k].mass;
			for (int j = 0; j < 3; j++)
			{
				massInvTriplet.push_back(Eigen::Triplet<double>(3 * (offset + k) + j, 3 * (offset + k) + j, massInverse));
			}
		}
		offset += nParticles;
	}


	// TODO: assemble rigid body mass inverse

	
	mInv.resize(3 * offset,3 * offset);
	mInv.setFromTriplets(massInvTriplet.begin(), massInvTriplet.end());
}



bool ElasticHook::numericalIntegration()
{
	// Unconstrained update
	int nbodies = (int)bodies_.size();

	std::vector<Eigen::Vector3d> oldthetas;
	for (int bodyidx = 0; bodyidx < (int)bodies_.size(); bodyidx++)
	{
		RigidBodyInstance &body = *bodies_[bodyidx];
		body.c += params_.timeStep*body.cvel;
		Eigen::Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*body.w);
		Eigen::Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

		Eigen::Vector3d oldtheta = body.theta;
		body.theta = VectorMath::axisAngle(Rtheta*Rhw);
		if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI / 2.0)
		{
			double oldnorm = oldtheta.norm();
			oldtheta = (oldnorm - 2.0*M_PI)*oldtheta / oldnorm;
		}
		body.oldtheta = oldtheta;
		oldthetas.push_back(oldtheta);
	}
	
	// Lagrangian, solved via Newton's method
	// TODO : Add Rigid body coupling constraint
	
	Eigen::VectorXd pos, vel;
	buildConfiguration(pos, vel);
	pos += params_.timeStep * vel;
	unbuildConfiguration(pos, vel);

	Eigen::VectorXd F;
	F.resizeLike(pos);
	F.setZero();
	computeForce(F);

	//std::cout << F.size() << std::endl;

	Eigen::SparseMatrix<double> gradG;
	Eigen::VectorXd g, lambda;
	computeContraintsAndGradient(pos, g, gradG);

	lambda.resizeLike(g);
	lambda.setZero();


	bool flag = newtonSolver(lambda, [this, &pos, &vel, &F](Eigen::VectorXd lambda , Eigen::VectorXd &f, Eigen::SparseMatrix<double> &gradF)
	{
		this->computeLagrangeMultiple(lambda, f, gradF, pos, vel, F);
	});
	if (flag)
	{
		// TODO: update rigid body
		std::cout << "lambda : " << lambda.norm() << std::endl;
		std::cout << "F : " << F.norm() << std::endl;
		std::cout << "gradG : " << gradG.norm() << std::endl;
		vel += (params_.timeStep*mInv_*(F + gradG.transpose()*lambda));
	}
	unbuildConfiguration(pos, vel);
	

	Eigen::VectorXd cForce(3 * nbodies);
	Eigen::VectorXd thetaForce(3 * nbodies);
	cForce.setZero();
	thetaForce.setZero();
	//computeForces(cForce, thetaForce);

	for (int bodyidx = 0; bodyidx < (int)bodies_.size(); bodyidx++)
	{
		RigidBodyInstance &body = *bodies_[bodyidx];
		Eigen::Matrix3d Mi = body.getTemplate().getInertiaTensor();

		body.cvel += params_.timeStep*cForce.segment<3>(3 * bodyidx) / body.density / body.getTemplate().getVolume();

		Eigen::Vector3d newwguess(body.w);

		int iter = 0;
		for (iter = 0; iter < params_.NewtonMaxIters; iter++)
		{
			Eigen::Vector3d term1 = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.oldtheta)).transpose() * Mi * body.density * newwguess;
			Eigen::Vector3d term2 = (VectorMath::TMatrix(params_.timeStep*body.w).inverse()*VectorMath::TMatrix(body.oldtheta)).transpose() * Mi * body.density * body.w;
			Eigen::Vector3d term3 = -params_.timeStep * thetaForce.segment<3>(3 * bodyidx);
			Eigen::Vector3d fval = term1 + term2 + term3;
			if (fval.norm() / body.density / Mi.trace() <= params_.NewtonTolerance)
				break;

			Eigen::Matrix3d Df = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

			Eigen::Vector3d deltaw = Df.inverse() * (-fval);
			newwguess += deltaw;
		}
		//std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
		body.w = newwguess;

	}




	return flag;
}

void ElasticHook::computeForce(Eigen::VectorXd & F)
{
	
	int offset = 0;

	for (int i = 0; i < (int)rods_.size(); i++)
	{
		int nvars = rods_[i]->nodes.size() * 3;
		Eigen::VectorXd f;
		rods_[i]->computeCenterlineForces(f);
		for (int k = 0; k < f.size(); k++)
		{
			F(offset + k) += f(k);
		}
		offset += nvars;
	} 

	// Compute External force here:
	if (params_.gravityEnabled)
		computeGravityForce(F);

	//TODO: Assemble Rigid body forces

}



void ElasticHook::computeGravityForce(Eigen::VectorXd & F)
{
	int offset = 0;
	Eigen::Vector3d gravityDir(0, 1.0, 0);
	for (int i = 0; i < (int)rods_.size(); i++)
	{
		int nParticles = rods_[i]->nodes.size();
		for (int k = 0; k < nParticles; k++)
		{
			F.segment<3>(3 * (offset + k)) += params_.gravityG * rods_[i]->nodes[k].mass * gravityDir;
		}
		offset += nParticles;
	}
}


void ElasticHook::computeContraintsAndGradient(const Eigen::VectorXd &q, Eigen::VectorXd &g, Eigen::SparseMatrix<double>& gradG)
{
	int nConfigurations = 0, nConstraints = 0;

	for (int i = 0; i < (int)rods_.size(); i++)
	{
		nConfigurations += (3 * rods_[i]->nodes.size());
		nConstraints += rods_[i]->nodes.size()-1;
	}

	g.resize(nConstraints);


	std::vector<Eigen::Triplet<double> > dgTriplet;

	// Inextensible Constarint
	int particleOffset = 0;
	int rodsOffset = 0;
	for (int rodidx= 0; rodidx < (int)rods_.size(); rodidx++)
	{
		ElasticRod &rod = *rods_[rodidx];
		int nparticles = (int)rod.nodes.size();
		for (int i = 0; i < nparticles - 1; i++)
		{
			Eigen::Vector3d e1, e2;
			e1 = q.segment<3>(3 * (particleOffset + i));
			e2 = q.segment<3>(3 * (particleOffset + i + 1));
			g(rodsOffset + i) = (e2 - e1).squaredNorm() - rod.restLength(i) * rod.restLength(i);

			Eigen::Vector3d localF = 2 * (e1 - e2);
			for (int k = 0; k < 3; k++)
			{
				dgTriplet.push_back(Eigen::Triplet<double>(rodsOffset + i, 3 * (particleOffset + i) + k, localF(k)));
				dgTriplet.push_back(Eigen::Triplet<double>(rodsOffset + i, 3 * (particleOffset + i + 1) + k, -localF(k)));

			}
		}
		particleOffset += nparticles;
		rodsOffset += (nparticles - 1);
	}
	
	// TODO: deal with rigid body coupling constarint

	gradG.resize(nConstraints, nConfigurations);
	gradG.setFromTriplets(dgTriplet.begin(), dgTriplet.end());
}

void ElasticHook::computeLagrangeMultiple(Eigen::VectorXd &lambda, Eigen::VectorXd &f, Eigen::SparseMatrix<double>&gradF, const Eigen::VectorXd &pos, const Eigen::VectorXd &vel, const Eigen::VectorXd &constF)
{

	Eigen::VectorXd g;
	Eigen::SparseMatrix<double> gradG;
	computeContraintsAndGradient(pos, g, gradG);
	
	// Change pos
	Eigen::VectorXd nextPos, nextVel;
	//std::cout << mInv_.size() << ' ' << constF.size() << std::endl;
	nextVel = vel + params_.timeStep*mInv_*(constF + gradG.transpose()*lambda);
	nextPos = pos + params_.timeStep*nextVel;

	computeContraintsAndGradient(nextPos, f, gradF);
	Eigen::SparseMatrix<double> W;
	W = params_.timeStep * params_.timeStep * mInv_ *gradG.transpose();
	gradF = gradF * W;
}

bool ElasticHook::newtonSolver(Eigen::VectorXd &x, std::function<void(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> _computeFAndGradF)
{
	double ratio = 1e-6;
	for (int i = 0; i<params_.NewtonMaxIters; i++)
	{
		Eigen::VectorXd F, dx;
		Eigen::SparseMatrix<double> gradF;
		_computeFAndGradF(x, F, gradF);
		if (F.norm()<params_.NewtonTolerance)
		{
			std::cout<<"Lagrangian Newton's Iteration Converges in " << i << "steps!" <<std::endl;
			return true;
		}
		Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
		solver.compute(gradF);
		dx = solver.solve(-F);
		x = x + dx;
	}
	//std::cout<<"Maximun iteration reached !!"<<std::endl;
	return true; 
}




//////////////////////////////////////////////////////////////////////////////////////
////    Helper Function
/////////////////////////////////////////////////////////////////////////////////////

void ElasticHook::testProcess()
{
	/*
	Here will put whatever the test pipeline is
	*/

	if (rods_.size() == 0) return;
	
	//Test energy differential
	
	ElasticRod &rod = *rods_[0];
	ElasticRod *temprod = rods_[0];
	
	Eigen::VectorXd pos, vel, thetas;
	temprod->buildConfiguration(pos, vel);

	int nrods = (int)temprod->rods.size();
	thetas.resize(nrods);
	for (int i = 0; i < nrods; i++)
	{
		thetas[i] = temprod->getRodSegment(i).theta;
	}

	std::cout << "\n\nTest centerline force:\n";
	auto energyFunc = [temprod](const Eigen::VectorXd &q) {
		Eigen::VectorXd dummy,v;
		temprod->buildConfiguration(dummy, v);
		temprod->unbuildConfiguration(q, v);
		double e=temprod->computeTotalEnergy();
		temprod->unbuildConfiguration(dummy, v);
		return e;
	};

	auto forceFunc = [temprod](const Eigen::VectorXd &q, Eigen::VectorXd &f) {
		Eigen::VectorXd dummy, v;
		//std::cout << q.transpose() << std::endl;
		temprod->buildConfiguration(dummy, v);
		temprod->unbuildConfiguration(q, v);
		//std::cout << "Here!" << std::endl;
		temprod->computeCenterlineForces(f);
		temprod->unbuildConfiguration(dummy, v);
	};

	TestModule::testEnergyDifferential(pos, energyFunc, forceFunc);

	std::cout << "\n\nTest Theta Force\n";

	auto energyThetaFunc = [temprod](const Eigen::VectorXd &q) {
		Eigen::VectorXd tmp;
		tmp.resizeLike(q);
		for (int k = 0; k < (int)temprod->rods.size(); k++)
		{
			tmp[k] = temprod->getRodSegment(k).theta;
			temprod->rods[k].theta = q(k);
		}
		temprod->updateMaterialCurvature();
		double e=temprod->computeTotalEnergy();
		for (int k = 0; k < (int)temprod->rods.size(); k++)
		{
			temprod->rods[k].theta = tmp(k);
		}
		temprod->updateMaterialCurvature();
		return -e;
	};

	auto thetaForceFunc = [temprod](const Eigen::VectorXd &q, Eigen::VectorXd &f) {
		Eigen::VectorXd dummy1, dummy2, dummy3;
		temprod->computeEnergyThetaDifferentialAndHessian(q, f, dummy1, dummy2, dummy3);
	};

	TestModule::testEnergyDifferential(thetas, energyThetaFunc, thetaForceFunc);

	std::cout << "\n\nTest Theta Hessian\n";
	auto thetaFunc = [temprod](const Eigen::VectorXd &q, Eigen::VectorXd &f, Eigen::SparseMatrix<double> &H) {
		Eigen::VectorXd lowerH, centerH, upperH;
		temprod->computeEnergyThetaDifferentialAndHessian(q, f, lowerH, centerH, upperH);
		H.resize(q.size(), q.size());
		std::vector<Eigen::Triplet<double>> hTriplet;
		for (int i = 0; i < centerH.size(); i++)
		{
			if (i > 0) hTriplet.push_back(Eigen::Triplet<double>(i, i - 1, lowerH[i]));
			hTriplet.push_back(Eigen::Triplet<double>(i, i, centerH[i]));
			if (i < centerH.size() - 1) hTriplet.push_back(Eigen::Triplet<double>(i, i + 1, upperH[i]));
		}
		H.setFromTriplets(hTriplet.begin(), hTriplet.end());
	};

	TestModule::testEnergyHessian(thetas, thetaFunc);
	temprod = NULL;
}

void ElasticHook::loadMesh()
{
	/*
	Load Mesh for Rendering Purposes
	TODO: Load Ridid bodies
	*/
	for (RigidBodyInstance *rbi : bodies_)
		delete rbi;
	for (RigidBodyTemplate *rbt : templates_)
		delete rbt;
	bodies_.clear();
	bodies_.shrink_to_fit();
	templates_.clear();
	templates_.shrink_to_fit();

	std::string prefix;
	std::string meshFilename = std::string("meshes/sphere.obj");
	std::ifstream ifs(meshFilename);
	if (!ifs)
	{
		// run from the build directory?
		prefix = "../";
		meshFilename = prefix + meshFilename;
		ifs.open(meshFilename);
		if (!ifs)
		{
			prefix = "../";
			meshFilename = prefix + meshFilename;
			ifs.open(meshFilename);
			prefix = "../../";
			if (!ifs)
				return;
		}
	}
	if (ballTemplate_) delete ballTemplate_;
	ballTemplate_ = new RigidBodyTemplate(meshFilename, 0.01);
	meshFilename = prefix + std::string("meshes/box.obj");
	if(rodTemplate_) delete rodTemplate_;
	rodTemplate_ = new RigidBodyTemplate(meshFilename, 0.01);
	std::cout << "Initialize Complete!\n";
}

void ElasticHook::saveConfiguration(std::string filePath)
{
	//TODO: Save scenes I created
}

void ElasticHook::loadConfiguration(std::string filePath)
{
	//TODO: load scenes I created
}
