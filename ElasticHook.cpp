
#include "ElasticHook.h"
#include "TestModule.h"

#include <Eigen/SparseQR>



void ElasticHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen))
	{
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Load", ImVec2((w - p) / 2.f, 0)))
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
		}
	}
	if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Combo("Click Adds", (int *)&params_.clickMode, "Particles\0Toggle Rotation\0\0");
		ImGui::Combo("Connector Type", (int *)&params_.connectorType, "Flexible Rods\0\0");
	}
	if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
	{
		//ImGui::Combo("Constraint Handling", (int *)&params_.constraintHandling, "Penalty Method\0Step and Project\0Lagrange Multipliers\0\0");
		//ImGui::Combo("Integrator", (int *)&params_.integrator, "Velocity Verlet\0Implicit Midpoint\0\0");
		ImGui::InputDouble("Timestep", &params_.timeStep);
		ImGui::InputDouble("Newton Tolerance", &params_.NewtonTolerance);
		ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
		ImGui::InputDouble("Penalty Stiffness", &params_.penaltyStiffness);
	}
	if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
		ImGui::InputDouble("  Gravity g", &params_.gravityG);
		ImGui::Checkbox("Floor Enabled", &params_.floorEnabled);
		ImGui::Checkbox("Bending Enabled", &params_.bendingEnabled);
		ImGui::Checkbox("Twist Enabled", &params_.twistEnabled);
		//viewer.imgui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
	}


	if (ImGui::CollapsingHeader("New Particles", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Checkbox("Is Fixed", &params_.particleFixed);
		ImGui::InputDouble("Mass", &params_.particleMass);
	}

	if (ImGui::CollapsingHeader("New Rods", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::InputInt("Num Segments", &params_.rodSegments);
		ImGui::InputDouble("Density", &params_.rodDensity);
		ImGui::InputDouble("Stretching Stiffness", &params_.rodStretchingStiffness);
		ImGui::InputDouble("Bending Stiffness", &params_.rodBendingStiffness);
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

	// this is terrible. But, easiest to get up and running

	std::vector<Eigen::Vector3d> verts;
	std::vector<Eigen::Vector3d> vertexColors;
	std::vector<Eigen::Vector3i> faces;

	int idx = 0;

	double eps = 1e-4;


	if (params_.floorEnabled)
	{
		/*for (int i = 0; i < 5; i++)
		{
		vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
		}*/

		verts.push_back(Eigen::Vector3d(0, -1, 0));
		verts.push_back(Eigen::Vector3d(1, -1, 1));
		verts.push_back(Eigen::Vector3d(-1, -1, 1));
		verts.push_back(Eigen::Vector3d(-1, -1, -1));
		verts.push_back(Eigen::Vector3d(1, -1, -1));

		faces.push_back(Eigen::Vector3i(idx + 0, idx + 2, idx + 1));
		faces.push_back(Eigen::Vector3i(idx + 0, idx + 3, idx + 2));
		faces.push_back(Eigen::Vector3i(idx + 0, idx + 4, idx + 3));
		faces.push_back(Eigen::Vector3i(idx + 0, idx + 1, idx + 4));

		idx += 5;
	}


	double basescale = 0.1;
	//Push a test bar 
	int nverts = rodV.rows();
	int nfaces = rodF.rows();
	for (int i = 0; i < nverts; i++)
	{
		verts.push_back(basescale * rodV.row(i));

	}
	for (int i = 0; i < nfaces; i++)
	{
		faces.push_back(rodF.row(i).transpose() + idx * Eigen::Vector3i::Ones());
	}
	idx += nverts;

	/*int nparticles = particles_.size();

	for (int i = 0; i<nparticles; i++)
	{
		int nverts = ballV.rows();
		int nfaces = ballF.rows();
		Eigen::RowVector3d c = Eigen::RowVector3d(particles_[i].pos[0], particles_[i].pos[1], 0);
		for (int i = 0; i < nverts; i++)
		{
			verts.push_back(c + (ballV.row(i) - c)*baseradius);
		}
		for (int i = 0; i < nfaces; i++)
		{
			faces.push_back(ballF.row(i).transpose() + idx * Eigen::Vector3i::Ones());
		}

		idx += nverts;

	}*/

	renderQ.resize(verts.size(), 3);
	//renderC.resize(vertexColors.size(), 3);
	for (int i = 0; i < verts.size(); i++)
	{
		renderQ.row(i) = verts[i];
		//renderC.row(i) = vertexColors[i];
	}
	renderF.resize(faces.size(), 3);
	for (int i = 0; i < faces.size(); i++)
		renderF.row(i) = faces[i];
}


void ElasticHook::initSimulation()
{
	time_ = 0;
	//particles_.clear();
	for (std::vector<ElasticRod *>::iterator it = rods_.begin(); it != rods_.end(); ++it)
		delete *it;
	loadMesh();
}

bool ElasticHook::mouseClicked(igl::opengl::glfw::Viewer & viewer, Eigen::Vector3d dir, int button)
{
	if (button != 2)
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
		while (!mouseClicks_.empty())
		{
			MouseClick mc = mouseClicks_.front();
			mouseClicks_.pop_front();
			switch (mc.mode)
			{
			case SimParameters::ClickMode::CM_ADDPARTICLE:
			{
				addParticle(mc.x, mc.y, 0);
				break;
			}
			}
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
	bool isSucceed = numericalIntegration();
	time_ += params_.timeStep;
	return (false || !isSucceed);
}

void ElasticHook::addParticle(double x, double y, double z)
{
	Eigen::Vector3d newpos(x, y, z);
	double mass = params_.particleMass;
	if (params_.particleFixed)
		mass = std::numeric_limits<double>::infinity();
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
		oldthetas.push_back(oldtheta);
	}

	int nrods = (int)rods_.size();
	for (int k = 0; k < nrods; k++)
	{
		ElasticRod &rod = *rods_[k];
		Eigen::VectorXd pos, vel;
		rod.buildConfiguration(pos, vel);
		pos += params_.timeStep * vel;
		rod.unbuildConfiguration(pos, vel);	
	}


	
	// Lagrangian, solved via Newton's method
	
	int nrods = (int)rods_.size();
	for (int k = 0; k < nrods; k++)
	{
		ElasticRod &rod = *rods_[k];
		
		Eigen::SparseMatrix<double> gradG;
		Eigen::VectorXd g, lambda;
		rod.computeInextensibleConstraint(g, gradG);

		Eigen::VectorXd F;
		rod.computeCenterlineForces(F);

		Eigen::VectorXd lambda;
		int nRods = (int)rod.nodes.size() - 1;
		lambda.resize(nRods);
		lambda.setZero();

		bool flag = newtonSolver(lambda, rod.computeLagrangeMultiple);
		if (flag)
		{
			Eigen::VectorXd vel = params_.timeStep * (F + gradG.transpose()*lambda);
			for (int i = 0; i < (int)(rod.nodes).size(); i++)
			{
				rod.nodes[i].vel += vel.segment<3>(3 * i) / rod.nodes[i].mass;
			}
		}
		else
		{
			std::cout << "Lagrangian blow up!\n" << std::endl;
		}

		rod.updateAfterTimeIntegration();

	}

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
			Eigen::Vector3d term1 = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
			Eigen::Vector3d term2 = (VectorMath::TMatrix(params_.timeStep*body.w).inverse()*VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
			Eigen::Vector3d term3 = -params_.timeStep * thetaForce.segment<3>(3 * bodyidx);
			Eigen::Vector3d fval = term1 + term2 + term3;
			if (fval.norm() / body.density / Mi.trace() <= params_.NewtonTolerance)
				break;

			Eigen::Matrix3d Df = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

			Eigen::Vector3d deltaw = Df.inverse() * (-fval);
			newwguess += deltaw;
		}
		//        std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
		body.w = newwguess;

	}




	return false;
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
			//std::cout<<"Optimal station reached!!"<<std::endl;
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


void ElasticHook::computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H)
{
	
}



//////////////////////////////////////////////////////////////////////////////////////
////    Helper Function
/////////////////////////////////////////////////////////////////////////////////////

void ElasticHook::testForceDifferential()
{

	/*Eigen::VectorXd q, vel, qPrev;
	buildConfiguration(q, vel, qPrev);

	// Note: Now previous position will not add to buildConfiguration since we do not need it.
	// If we have to implement some other numerical integrator we need to change buildConfiguration and
	// unbuildConfiguration
	qPrev.resize(q.size());
	for (int i = 0; i<(int)particles_.size(); i++)
	{
	qPrev.segment<2>(2 * i) = particles_[i].prevpos;
	}

	Eigen::SparseMatrix<double> gradF(q.size(), q.size());
	gradF.setZero();
	Eigen::VectorXd f = Eigen::VectorXd::Zero(q.size());
	computeForceAndHessian(q, qPrev, f, gradF);

	Eigen::VectorXd direction = Eigen::VectorXd::Random(q.size());
	for(int i=0;i<particles_.size();i++)
	{
	if(particles_[i].fixed)
	direction.segment(2*i, 2).setZero();
	}
	direction.normalize();

	for (int k = 1; k <= 12; k++)
	{
	double eps = pow(10, -k);

	VectorXd epsF = Eigen::VectorXd::Zero(q.size());
	Eigen::SparseMatrix<double>  epsgradF(q.size(), q.size());
	computeForceAndHessian(q + eps * direction, qPrev, epsF, epsgradF);

	std::cout << "EPS is: " << eps << std::endl;
	std::cout << "Norm of Finite Difference is: " << (epsF - f).norm() / eps << std::endl;
	std::cout << "Norm of Directinal Gradient is: " << (gradF*direction).norm() << std::endl;
	std::cout << "The difference between above two is: " << ((epsF - f) / eps + gradF*direction).norm() << std::endl << std::endl;
	}*/
}




void ElasticHook::loadMesh()
{
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

	igl::readOBJ(meshFilename, ballV, ballF);
	meshFilename = prefix + std::string("meshes/box.obj");
	igl::readOBJ(meshFilename, rodV, rodF);
	std::cout << "Initialize Complete!\n";
}

void ElasticHook::saveConfiguration(std::string filePath)
{
	/*int nParticles = particles_.size();
	int nConnectors = connectors_.size();
	int nBendingStencils = bendingStencils_.size();
	std::ofstream outfile(filePath, std::ios::trunc);

	outfile << nParticles << "\n";
	outfile << nConnectors << "\n";
	outfile << nBendingStencils << "\n";

	for (int i = 0; i<nParticles; i++)
	{
	outfile << std::setprecision(16) << particles_[i].pos(0) << "\n";
	outfile << std::setprecision(16) << particles_[i].pos(1) << "\n";
	outfile << std::setprecision(16) << particles_[i].prevpos(0) << "\n";
	outfile << std::setprecision(16) << particles_[i].prevpos(1) << "\n";
	outfile << std::setprecision(16) << particles_[i].vel(0) << "\n";
	outfile << std::setprecision(16) << particles_[i].vel(1) << "\n";

	outfile << particles_[i].fixed << "\n";
	if (!particles_[i].fixed)
	outfile << std::setprecision(16) << particles_[i].mass << "\n";
	outfile << particles_[i].inert << "\n";
	}

	for (int i = 0; i<nConnectors; i++)
	{
	outfile << connectors_[i]->getType()<<"\n";
	outfile << connectors_[i]->p1 << "\n";
	outfile << connectors_[i]->p2 << "\n";
	outfile << connectors_[i]->mass << "\n";
	switch (connectors_[i]->getType())
	{
	case SimParameters::CT_SPRING:
	{
	Spring &s = *(Spring *) connectors_[i];
	outfile << s.canSnap << "\n";
	outfile << std::setprecision(16) << s.stiffness << "\n";
	break;
	}
	case SimParameters::CT_RIGIDROD:
	{
	break;
	}
	default:
	std::cout << "Unknown Connector type, your saving might fail.\n";
	break;
	}
	outfile << connectors_[i]->associatedBendingStencils.size() << "\n";

	for (std::set<int>::iterator sit = connectors_[i]->associatedBendingStencils.begin(); sit != connectors_[i]->associatedBendingStencils.end(); ++sit)
	{
	outfile << (*sit) << "\n";
	}

	}

	for (int i = 0; i < nBendingStencils; i++)
	{
	outfile << bendingStencils_[i].p1 << "\n";
	outfile << bendingStencils_[i].p2 << "\n";
	outfile << bendingStencils_[i].p3 << "\n";
	outfile << std::setprecision(16) << bendingStencils_[i].kb << "\n";
	//outfile << std::setprecision(16) << bendingStencils_[i].theta << "\n";
	}
	outfile.close();*/
}

void ElasticHook::loadConfiguration(std::string filePath)
{
	/*std::ifstream infile(filePath);
	if (!infile)
	return;
	int nParticles;
	int nConnectors;
	int nBendingStencils;

	infile >> nParticles;
	infile >> nConnectors;
	infile >> nBendingStencils;

	particles_.clear();
	// Only available after C++ 11, hopefully won't throw compile error in GDC's computer
	particles_.shrink_to_fit();
	connectors_.clear();
	connectors_.shrink_to_fit();
	bendingStencils_.clear();
	bendingStencils_.shrink_to_fit();


	for (int i = 0; i < nParticles; i++)
	{
	Eigen::Vector2d pos, prevpos, vel;
	infile >> pos(0);
	infile >> pos(1);
	infile >> prevpos(0);
	infile >> prevpos(1);
	infile >> vel(0);
	infile >> vel(1);
	double mass;
	bool isFixed;
	bool isInert;

	infile >> isFixed;
	if (isFixed)
	mass = mass = std::numeric_limits<double>::infinity();
	else
	infile >> mass;
	infile >> isInert;
	Particle newParticle = Particle(pos, mass, isFixed, isInert);
	newParticle.prevpos = prevpos;
	newParticle.vel = vel;
	particles_.push_back(newParticle);
	}

	for (int i = 0; i < nConnectors; i++)
	{
	int p1, p2, type;
	double mass;
	infile >> type;
	infile >> p1;
	infile >> p2;
	infile >> mass;
	double dist = (particles_[p1].pos - particles_[p2].pos).norm();

	Connector* connector;
	switch (type)
	{
	case SimParameters::CT_SPRING:
	{
	bool isSnappable;
	double stiffness;
	infile >> isSnappable;
	infile >> stiffness;
	connector = new Spring(p1, p2, mass, stiffness, dist, isSnappable);
	break;
	}
	case SimParameters::CT_RIGIDROD:
	{
	connector = new RigidRod(p1, p2, mass, dist);
	break;
	}
	default:
	std::cout << "Unknown Connector type, your loading might fail.\n";
	break;
	}
	int nAssociate;
	infile >> nAssociate;
	for (int j = 0; j < nAssociate; j++)
	{
	int id;
	infile >> id;
	connector->associatedBendingStencils.insert(id);
	}
	connectors_.push_back(connector);
	}

	for (int i = 0; i<nBendingStencils; i++)
	{
	int p1, p2, p3;
	double bendingStiffness, theta;
	infile >> p1;
	infile >> p2;
	infile >> p3;
	infile >> bendingStiffness;
	//infile >> theta;
	bendingStencils_.push_back(BendingStencil(p1, p2, p3, bendingStiffness));
	bendingStencils_[i].theta = theta;
	}*/
}
