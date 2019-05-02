#include "GooHook.h"
#include <Eigen/SparseQR>
using namespace Eigen;


void GooHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	if (ImGui::CollapsingHeader("Configuration", ImGuiTreeNodeFlags_DefaultOpen))
	{
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Load", ImVec2((w - p) / 2.f, 0)))
		{
			if(!isPaused())
				pause();
			//auto gooHook = static_cast<GooHook*>(hook);
			std::string filePath = igl::file_dialog_open();
			loadConfiguration(filePath);
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Save", ImVec2((w - p) / 2.f, 0)))
		{
			//auto gooHook = static_cast<GooHook*>(hook);
			std::string filePath = igl::file_dialog_save();
			saveConfiguration(filePath);
		}
	}
    if (ImGui::CollapsingHeader("UI Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Combo("Click Adds", (int *)&params_.clickMode, "Particles\0\0");
        ImGui::Combo("Connector Type", (int *)&params_.connectorType, "Flexible Rods\0\0");
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        //ImGui::Combo("Constraint Handling", (int *)&params_.constraintHandling, "Penalty Method\0Step and Project\0Lagrange Multipliers\0\0");
		//ImGui::Combo("Integrator", (int *)&params_.integrator, "Velocity Verlet\0Implicit Midpoint\0\0");
        ImGui::InputDouble("Timestep",  &params_.timeStep);
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


void GooHook::updateRenderGeometry()
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


    if(params_.floorEnabled)
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
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        switch((*it)->getType())
        {
        case SimParameters::CT_SPRING:
		case SimParameters::CT_RIGIDROD:
		case SimParameters::CT_ELASTICROD:
		{
			Vector3d sourcepos = particles_[(*it)->p1].pos;
			Vector3d destpos = particles_[(*it)->p2].pos;

			Vector3d vec = destpos - sourcepos;
			Vector3d width = Eigen::Vector3d(vec.norm(),0.08,0.10);
			Eigen::Vector3d c = Eigen::Vector3d((sourcepos[0] + destpos[0]) / 2.0, (sourcepos[1] + destpos[1]) / 2.0, 0);
			Eigen::Vector3d rot = Eigen::Vector3d(0, 0, std::atan2(vec[1],vec[0]));

			int nverts = rodV.rows();
			int nfaces = rodF.rows();
			for (int i = 0; i < nverts; i++)
			{
				Eigen::Vector3d pos = c + basescale * VectorMath::rotationMatrix(rot) * width.cwiseProduct(rodV.row(i).transpose());
				verts.push_back(pos);
			}
			for (int i = 0; i < nfaces; i++)
			{
				faces.push_back(rodF.row(i).transpose() +idx * Eigen::Vector3i::Ones());
			}

			idx += nverts;
            break;
        }
        default:
            break;
        }
    }

    int nparticles = particles_.size();

    for(int i=0; i<nparticles; i++)
    {
		int nverts = ballV.rows();
		int nfaces = ballF.rows();
		Eigen::RowVector3d c=Eigen::RowVector3d(particles_[i].pos[0], particles_[i].pos[1], 0);
		for (int i = 0; i < nverts; i++)
		{
			verts.push_back(c+ (ballV.row(i)-c)*baseradius);
		}
		for (int i = 0; i < nfaces; i++)
		{
			faces.push_back(ballF.row(i).transpose() + idx * Eigen::Vector3i::Ones());
		}

		idx += nverts;
        
    }

    

    renderQ.resize(verts.size(),3);
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


void GooHook::initSimulation()
{
    time_ = 0;
    particles_.clear();
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
        delete *it;
    connectors_.clear();
    bendingStencils_.clear();
	loadMesh();
}

void GooHook::tick()
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
                addParticle(mc.x, mc.y);
                break;
            }
            }
        }
    }
    message_mutex.unlock();
}

bool GooHook::simulateOneStep()
{
  
    VectorXd q, v, prevq;
    buildConfiguration(q, v, prevq);
    bool isSucceed = numericalIntegration(q, v, prevq);
    unbuildConfiguration(q, v);
    time_ += params_.timeStep;
    return (false || !isSucceed);
}

void GooHook::addParticle(double x, double y)
{
    Vector3d newpos(x,y,0);
    double mass = params_.particleMass;
    if(params_.particleFixed)
        mass = std::numeric_limits<double>::infinity();

    int newid = particles_.size();
    particles_.push_back(Particle(newpos, mass, params_.particleFixed, false));
	int numparticles = particles_.size()-1;

    for(int i=0; i<numparticles; i++)
    {
        if(particles_[i].inert)
            continue;
        Vector3d pos = particles_[i].pos;
        double dist = (pos-newpos).norm();
        if(dist <= params_.maxSpringDist)
        {
            switch(params_.connectorType)
            {
            case SimParameters::CT_RIGIDROD:
            {
                connectors_.push_back(new RigidRod(newid, i, 0, dist));
                break;
            }
            case SimParameters::CT_FLEXROD:
			{
				// 1. Add max{rodSegments, 2}-1 particles in between the rod
				// 2. Connect them consequtively as springs (new_id, current_id, current_id+1, ..., i)
				// 3. Add stencils for each pair of the spring
				int prevprevid = newid;
				int previd = newid;
				int s = (params_.rodSegments > 2 ? params_.rodSegments : 2);
				double rodrestlen = dist / s;
				/*
				for (int j = 1; j <= s; j++)
				{
					Eigen::Vector2d hingepos = j * 1.0 / s * pos + (s - j) * 1.0 / s * newpos;
					int currentid = (int)particles_.size();
					// Deal with last segment carefully!
					if (j < s)
					{
						particles_.push_back(Particle(hingepos, 0.0, false, true));
					}
					else
					{
						currentid = i;
					}
					// Add springs
					int nconnectors = (int)connectors_.size();
					connectors_.push_back(new Spring(previd, currentid, params_.rodDensity * rodrestlen
						, params_.rodStretchingStiffness / rodrestlen, rodrestlen, false));
					if (j > 1)
					{
						// Add stencils	
						Spring &s1 = *(Spring *)connectors_[nconnectors];
						Spring &s2 = *(Spring *)connectors_[nconnectors - 1];
						double rodrealstiffness = params_.rodBendingStiffness / (s1.restlen + s2.restlen);
						int bendidx = (int)bendingStencils_.size();
						bendingStencils_.push_back(BendingStencil(prevprevid, previd, currentid, rodrealstiffness));
						//Associate stencil with two springs
						connectors_[nconnectors]->associatedBendingStencils.insert(bendidx);
						connectors_[nconnectors - 1]->associatedBendingStencils.insert(bendidx);
					}
					prevprevid = previd;
					previd = currentid;
				}*/
                break;
            }
            default:
                break;
            }
        }
    }
}

double GooHook::getTotalParticleMass(int idx)
{
	return 0.0;
}

void GooHook::buildConfiguration(VectorXd &q, VectorXd &v, VectorXd &prevq)
{
   
}

void GooHook::unbuildConfiguration(const VectorXd &q, const VectorXd &v)
{
   
}

int GooHook::getNumRigidRods()
{
	return 0;
}

bool GooHook::numericalIntegration(VectorXd &q, VectorXd &v, VectorXd &prevq)
{
	return false;
}


bool GooHook::newtonSolver(Eigen::VectorXd &x, std::function<void (Eigen::VectorXd, Eigen::VectorXd &, Eigen::SparseMatrix<double> *)> _computeFAndGradF)
{
    double ratio = 1e-6;
    for(int i=0; i<params_.NewtonMaxIters; i++)
    {
//        bool isSuccess = takeOneStep(ratio, x, x0, _computeFAndGradF);
//        if(!isSuccess)
//        {
//            std::cout<<"Failed to update the position"<<std::endl;
//            return false;
//        }
//        Eigen::VectorXd F;
//        _computeFAndGradF(x, x0, F, NULL);
//        if(F.norm()<params_.NewtonTolerance)
//        {
//            std::cout<<"Optimal station reached!!"<<std::endl;
//            return true;
//        }
        Eigen::VectorXd F,dx;
        Eigen::SparseMatrix<double> gradF;
        _computeFAndGradF(x, F, &gradF);
        if(F.norm()<params_.NewtonTolerance)
        {
            //std::cout<<"Optimal station reached!!"<<std::endl;
            return true;
        }
        Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        solver.compute(gradF);
        dx = solver.solve(-F);
        x = x + dx;
    }
    //std::cout<<"Maximun iteration reached !!"<<std::endl;
    return true;
}

bool GooHook::takeOneStep(double &ratio, Eigen::VectorXd &x, std::function<void (Eigen::VectorXd,  Eigen::VectorXd &, Eigen::SparseMatrix<double> *)> _computeFAndGradF)
{
    int i;
    for(i=0;i<1e2;i++)
    {
        Eigen::VectorXd F,dx;
        Eigen::SparseMatrix<double> gradF, I;
        _computeFAndGradF(x, F, &gradF);
        I.resize(x.size(), x.size());
        I.setIdentity();
        gradF += ratio * I;
        if(F.norm()<params_.NewtonTolerance)
            break;
        Eigen::SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        solver.compute(gradF);
        dx = solver.solve(-F);
        Eigen::VectorXd newX = x + dx;
        
        Eigen::VectorXd newF;
        _computeFAndGradF(newX, newF, NULL);
        
        
        if(newF.norm() <= F.norm())
        {
            //std::cout<<"Sucessful update: "<<std::endl;
            //std::cout<<"Position Change: "<<dx.norm()<<" Norm of old force: "<<F.norm()<<" Norm of new force: "<<newF.norm()<<std::endl;
            x = newX;
            ratio = std::max(ratio / 2, 1e-6);
            break;
        }
        else
        {
            
            //std::cout<<"Norm of old force: "<<F.norm()<<" Norm of new force: "<<newF.norm()<<" lambda now: "<<ratio<<std::endl;
            ratio = ratio * 2;
        }
    }
    if(i<1e2)
        return true;
    else
        return false;
}


void GooHook::computeForceAndHessian(const VectorXd &q, const VectorXd &qprev, Eigen::VectorXd &F, SparseMatrix<double> &H)
{
}

void GooHook::processGravityForce(VectorXd &F)
{
}

void GooHook::processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
 
}

void GooHook::processDampingForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
 
}

void GooHook::processFloorForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
}

void GooHook::processPenaltyForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    
}

void GooHook::processBendingForce(const Eigen::VectorXd & q, const Eigen::VectorXd & qprev, Eigen::VectorXd & F, std::vector<Eigen::Triplet<double>>& H)
{
	int nBendingStencils = (int)bendingStencils_.size();
	for (int i = 0; i < nBendingStencils; i++)
	{
		BendingStencil hinge = bendingStencils_[i];
		int idx1 = hinge.p1;
		int idx2 = hinge.p2;
		int idx3 = hinge.p3;

		Vector2d p1 = q.segment<2>(2 * idx1);
		Vector2d p2 = q.segment<2>(2 * idx2);
		Vector2d p3 = q.segment<2>(2 * idx3);
		Vector2d direction1 = p2 - p1;
		Vector2d direction2 = p3 - p2;
		double dist1 = direction1.norm();
		double dist2 = direction2.norm();

		double theta = 2 * std::atan2(direction1[0] * direction2[1] - direction2[0] * direction1[1] , dist1 * dist2 + direction1.dot(direction2));
		bendingStencils_[i].theta = theta;
		double coef1 = hinge.kb * theta / dist1 / dist1;
		double coef2 = hinge.kb * theta / dist2 / dist2;

		Vector2d F1(0.0, 0.0);
		F1[0] = coef1 * direction1[1];
		F1[1] = -coef1 * direction1[0];

		Vector2d F2(0.0, 0.0);
		F2[0] = coef2 * direction2[1];
		F2[1] = -coef2 * direction2[0];

		F.segment<2>(2 * idx1) += F1;
		F.segment<2>(2 * idx2) -= F1;

		F.segment<2>(2 * idx3) += F2;
		F.segment<2>(2 * idx2) -= F2;


		// Hessian H = -dF
		double denominator = dist1 * dist2 + direction1.dot(direction2);
		double numerator = direction1[0] * direction2[1] - direction2[0] * direction1[1];
		double coefdtheta1 = hinge.kb * 2 / (denominator * denominator + numerator * numerator) / dist1 / dist1;
		double coefdtheta2 = hinge.kb * 2 / (denominator * denominator + numerator * numerator) / dist2 / dist2;
		Eigen::Matrix2d J;
		J << 0, 1, 
			-1, 0;


		// J.transpose()=-J
		Eigen::Vector2d dthetad1 = -denominator * J * direction2 - numerator * (-direction1 * dist2 / dist1 - direction2);
		Eigen::Vector2d dthetad3 = denominator * J.transpose() * direction1 - numerator * (direction2 * dist1 / dist2 + direction1);
		Eigen::Vector2d dthetad2 = -dthetad1 - dthetad3;


		// \frac{\partial F_i}{\partial p_i}
		Eigen::Matrix2d localH1 = coefdtheta1 * dthetad1 * direction1.transpose() * J.transpose() 
			+ coef1 * ( - J.transpose() + 2 * direction1 * direction1.transpose() * J.transpose() / dist1 / dist1);
		// \frac{\partial F_k}{\partial p_i}
		Eigen::Matrix2d localH2 = coefdtheta2 * dthetad1 * direction2.transpose() * J.transpose();

		for (int j = 0; j<2; j++)
			for (int k = 0; k<2; k++)
			{
				H.push_back(Eigen::Triplet<double>(2 * idx1 + j, 2 * idx1 + k, -localH1.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx1 + j, 2 * idx2 + k, localH1.coeff(j, k)+ localH2.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx1 + j, 2 * idx3 + k, -localH2.coeff(j, k)));
			}

		localH1 = coefdtheta1 * dthetad2 * direction1.transpose() * J.transpose()
			+ coef1 * (J.transpose() - 2 * direction1 * direction1.transpose() * J.transpose() / dist1 / dist1);

		localH2 = coefdtheta2 * dthetad2 * direction2.transpose() * J.transpose()
			+ coef2 * (-J.transpose() + 2 * direction2 * direction2.transpose() * J.transpose() / dist2 / dist2);

		for (int j = 0; j<2; j++)
			for (int k = 0; k<2; k++)
			{
				H.push_back(Eigen::Triplet<double>(2 * idx2 + j, 2 * idx1 + k, -localH1.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx2 + j, 2 * idx2 + k, localH1.coeff(j, k) + localH2.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx2 + j, 2 * idx3 + k, -localH2.coeff(j, k)));
			}

		localH1 = coefdtheta1 * dthetad3 * direction1.transpose() * J.transpose();
		localH2 = coefdtheta2 * dthetad3 * direction2.transpose() * J.transpose() 
			+ coef2 * (J.transpose() - 2 * direction2 * direction2.transpose() * J.transpose() / dist2 / dist2);

		for (int j = 0; j<2; j++)
			for (int k = 0; k<2; k++)
			{
				H.push_back(Eigen::Triplet<double>(2 * idx3 + j, 2 * idx1 + k, -localH1.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx3 + j, 2 * idx2 + k, localH1.coeff(j, k)+localH2.coeff(j, k)));
				H.push_back(Eigen::Triplet<double>(2 * idx3 + j, 2 * idx3 + k, -localH2.coeff(j, k)));
			}

	}
}





//////////////////////////////////////////////////////////////////////////////////////
////    Helper Function
/////////////////////////////////////////////////////////////////////////////////////

void GooHook::testForceDifferential()
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




void GooHook::loadMesh()
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

void GooHook::saveConfiguration(std::string filePath)
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

void GooHook::loadConfiguration(std::string filePath)
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
