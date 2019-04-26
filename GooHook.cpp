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
        ImGui::Combo("Click Adds", (int *)&params_.clickMode, "Particles\0Saws\0\0");
        ImGui::Combo("Connector Type", (int *)&params_.connectorType, "Springs\0Rigid Rods\0Flexible Rods\0\0");
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Combo("Constraint Handling", (int *)&params_.constraintHandling, "Penalty Method\0Step and Project\0Lagrange Multipliers\0\0");
		ImGui::Combo("Integrator", (int *)&params_.integrator, "Velocity Verlet\0Implicit Midpoint\0\0");
        ImGui::InputDouble("Timestep",  &params_.timeStep);
        ImGui::InputDouble("Newton Tolerance", &params_.NewtonTolerance);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
        ImGui::InputDouble("Penalty Stiffness", &params_.penaltyStiffness);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputDouble("  Gravity g", &params_.gravityG);
        ImGui::Checkbox("Springs Enabled", &params_.springsEnabled);
        ImGui::InputDouble("  Max Strain", &params_.maxSpringStrain);
        ImGui::Checkbox("Damping Enabled", &params_.dampingEnabled);
        ImGui::InputDouble("  Viscosity", &params_.dampingStiffness);
        ImGui::Checkbox("Floor Enabled", &params_.floorEnabled);
        ImGui::Checkbox("Bending Enabled", &params_.bendingEnabled);
        //viewer.imgui->addWindow(Eigen::Vector2i(1000, 0), "New Objects");
    }


    if (ImGui::CollapsingHeader("New Particles", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Is Fixed", &params_.particleFixed);
        ImGui::InputDouble("Mass", &params_.particleMass);
    }

    if (ImGui::CollapsingHeader("New Saws", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Radius", &params_.sawRadius);
    }

    if (ImGui::CollapsingHeader("New Springs", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputDouble("Max Spring Dist", &params_.maxSpringDist);
        ImGui::InputDouble("Base Stiffness", &params_.springStiffness);
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
    double baseradius = 0.02;
    double pulsefactor = 0.1;
    double pulsespeed = 50.0;

    int sawteeth = 20;
    double sawdepth = 0.1;
    double sawangspeed = 10.0;

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
        for (int i = 0; i < 6; i++)
        {
            vertexColors.push_back(Eigen::Vector3d(0.3, 1.0, 0.3));
        }

        verts.push_back(Eigen::Vector3d(-1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(-1, -1, eps));

        faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));

        verts.push_back(Eigen::Vector3d(-1, -1, eps));
        verts.push_back(Eigen::Vector3d(1, -0.5, eps));
        verts.push_back(Eigen::Vector3d(1, -1, eps));
        faces.push_back(Eigen::Vector3i(idx + 3, idx + 4, idx + 5));
        idx += 6;
    }


    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        switch((*it)->getType())
        {
        case SimParameters::CT_SPRING:
        {
            Eigen::Vector3d color;
            if((*it)->associatedBendingStencils.empty())
                color << 0.0, 0.0, 1.0;
            else
                color << 0.75, 0.5, 0.75;
            Vector2d sourcepos = particles_[(*it)->p1].pos;
            Vector2d destpos   = particles_[(*it)->p2].pos;

            Vector2d vec = destpos - sourcepos;
            Vector2d perp(-vec[1], vec[0]);
            perp /= perp.norm();

            double dist = (sourcepos-destpos).norm();

            double width = baselinewidth/(1.0+ 20.0 * dist * dist);

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;

            break;
        }
        case SimParameters::CT_RIGIDROD:
        {
            Eigen::Vector3d color;
            if((*it)->associatedBendingStencils.empty())
                color << 1.0, 0.0, 1.0;
            else
                color << 1.0, 1.0, 0.3;

            Vector2d sourcepos = particles_[(*it)->p1].pos;
            Vector2d destpos   = particles_[(*it)->p2].pos;
            Vector2d vec = destpos - sourcepos;
            Vector2d perp(-vec[1], vec[0]);
            perp /= perp.norm();

            double width = baselinewidth;

            for (int i = 0; i < 4; i++)
                vertexColors.push_back(color);

            verts.push_back(Eigen::Vector3d(sourcepos[0] + width * perp[0], sourcepos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(sourcepos[0] - width * perp[0], sourcepos[1] - width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] + width * perp[0], destpos[1] + width * perp[1], -eps));
            verts.push_back(Eigen::Vector3d(destpos[0] - width * perp[0], destpos[1] - width * perp[1], -eps));

            faces.push_back(Eigen::Vector3i(idx, idx + 1, idx + 2));
            faces.push_back(Eigen::Vector3i(idx + 2, idx + 1, idx + 3));
            idx += 4;

            break;
        }
        default:
            break;
        }
    }

    int nparticles = particles_.size();

    for(int i=0; i<nparticles; i++)
    {
        double radius = baseradius*sqrt(getTotalParticleMass(i));
        radius *= (1.0 + pulsefactor*sin(pulsespeed*time_));

        Eigen::Vector3d color(0,0,0);

        if(particles_[i].fixed)
        {
            radius = baseradius;
            color << 1.0, 0, 0;
        }

        for (int j = 0; j < numcirclewedges + 2; j++)
        {
            vertexColors.push_back(color);
        }


        verts.push_back(Eigen::Vector3d(particles_[i].pos[0], particles_[i].pos[1], 0));

        const double PI = 3.1415926535898;
        for (int j = 0; j <= numcirclewedges; j++)
        {
            verts.push_back(Eigen::Vector3d(particles_[i].pos[0] + radius * cos(2 * PI*j / numcirclewedges),
                particles_[i].pos[1] + radius * sin(2 * PI*j / numcirclewedges), 0));
        }

        for (int j = 0; j <= numcirclewedges; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (numcirclewedges + 1))));
        }

        idx += numcirclewedges + 2;
    }

    for(std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
    {
        double outerradius = it->radius;
        double innerradius = (1.0-sawdepth)*outerradius;

        Eigen::Vector3d color(0.5,0.5,0.5);

        int spokes = 2*sawteeth;
        for (int j = 0; j < spokes + 2; j++)
        {
            vertexColors.push_back(color);
        }

        verts.push_back(Eigen::Vector3d(it->pos[0], it->pos[1], 0));

        const double PI = 3.1415926535898;
        for (int i = 0; i <= spokes; i++)
        {
            double radius = (i % 2 == 0) ? innerradius : outerradius;
            verts.push_back(Eigen::Vector3d(it->pos[0] + radius * cos(2 * PI*i / spokes + sawangspeed*time_),
                it->pos[1] + radius * sin(2 * PI*i / spokes + sawangspeed*time_), 0));
        }

        for (int j = 0; j <= spokes; j++)
        {
            faces.push_back(Eigen::Vector3i(idx, idx + j + 1, idx + 1 + ((j + 1) % (spokes + 1))));
        }

        idx += spokes + 2;
    }

    renderQ.resize(verts.size(),3);
    renderC.resize(vertexColors.size(), 3);
    for (int i = 0; i < verts.size(); i++)
    {
        renderQ.row(i) = verts[i];
        renderC.row(i) = vertexColors[i];
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
    saws_.clear();
    bendingStencils_.clear();
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
            case SimParameters::ClickMode::CM_ADDSAW:
            {
                addSaw(mc.x, mc.y);
                break;
            }
            }
        }
    }
    message_mutex.unlock();
}

bool GooHook::simulateOneStep()
{
    // TODO handle constraints and flexible rods
    VectorXd q, v, prevq;
    buildConfiguration(q, v, prevq);
    bool isSucceed = numericalIntegration(q, v, prevq);
    unbuildConfiguration(q, v);

    pruneOverstrainedSprings();
    deleteSawedObjects();
    time_ += params_.timeStep;
    return (false || !isSucceed);
}

void GooHook::addParticle(double x, double y)
{
    Vector2d newpos(x,y);
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
        Vector2d pos = particles_[i].pos;
        double dist = (pos-newpos).norm();
        if(dist <= params_.maxSpringDist)
        {
            switch(params_.connectorType)
            {
            case SimParameters::CT_SPRING:
            {
                connectors_.push_back(new Spring(newid, i, 0, params_.springStiffness/dist, dist, true));
                break;
            }
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
				}
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
    double mass = particles_[idx].mass;
	// Fixed points are unaffected.
	if (particles_[idx].fixed)
		return mass;
	//Iterate over all possible connectors
	for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
	{
		if (idx == (*it)->p1 || idx == (*it)->p2)
		{
			mass = mass + (*it)->mass / 2;
		}
	}
    return mass;
}

void GooHook::addSaw(double x, double y)
{
    saws_.push_back(Saw(Vector2d(x,y), params_.sawRadius));
}

void GooHook::buildConfiguration(VectorXd &q, VectorXd &v, VectorXd &prevq)
{
    int ndofs = 2*particles_.size();
    q.resize(ndofs);
    v.resize(ndofs);
	prevq.resize(ndofs);

    for(int i=0; i<(int)particles_.size(); i++)
    {
        q.segment<2>(2*i) = particles_[i].pos;
        v.segment<2>(2*i) = particles_[i].vel;
		prevq.segment<2>(2 * i) = particles_[i].prevpos;
    }
    
    computeMassInverse(Minv_);
}

void GooHook::unbuildConfiguration(const VectorXd &q, const VectorXd &v)
{
    int ndofs = q.size();
    assert(ndofs == int(2*particles_.size()));

    for(int i=0; i<ndofs/2; i++)
    {
		// Add this line back
		particles_[i].prevpos = particles_[i].pos;
        particles_[i].pos = q.segment<2>(2*i);
        particles_[i].vel = v.segment<2>(2*i);
    }
}

int GooHook::getNumRigidRods()
{
    int nrods = 0;
    for(std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
    {
        if((*it)->getType() == SimParameters::CT_RIGIDROD)
            nrods++;
    }
    return nrods;
}

bool GooHook::numericalIntegration(VectorXd &q, VectorXd &v, VectorXd &prevq)
{
    // TODO handle constraints and flexible rods
    VectorXd F;
    SparseMatrix<double> H;

    computeMassInverse(Minv_);

    VectorXd oldq = q;
    VectorXd oldv = v;
    bool flag = true;
    if(params_.constraintHandling == params_.CH_PENALTY)
    {
		switch (params_.integrator)
		{
		case SimParameters::TI_VELOCITY_VERLET:
		{
			updatebyVelocityVerletUnconstraint(q,v);
			break;
		}
		case SimParameters::TI_IMPLICIT_MIDPOINT:
		{
			updatebyImpliciyMidpointUnconstraint(q, v, prevq);
			break;
		}
		default:
			break;
		}
    }
    else if(params_.constraintHandling == params_.CH_STEPPROJECT)
    {
        switch (params_.integrator)
        {
            case SimParameters::TI_VELOCITY_VERLET:
            {
                q += params_.timeStep*v;
                computeForceAndHessian(q, oldq, F, H);
                v += params_.timeStep*Minv_*F;
                break;
            }
            case SimParameters::TI_IMPLICIT_MIDPOINT:
            {
                
                // Update by midpoint
                flag = newtonSolver(q, [this, oldq, prevq, v](Eigen::VectorXd pos, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF)
                {
                    Eigen::SparseMatrix<double> idMat(pos.size(), pos.size());
                    idMat.setIdentity();
                    Eigen::VectorXd force;
                    Eigen::SparseMatrix<double> H;
                    
                    // H = -gradF
                    this->computeForceAndHessian((pos + oldq) / 2, (oldq + prevq) / 2, force, H);
                    
                    F = pos - oldq - params_.timeStep*v - 0.5*params_.timeStep*params_.timeStep*Minv_*force;
                    
                    *gradF = idMat - params_.timeStep * params_.timeStep * Minv_ * (-H) / 4.0;
                });
                if(flag)
                {
                    v = 2.0 * (q - oldq) / params_.timeStep - v;
                }
                else
                {
                    std::cout<<"Failed to update by implicit midpoint. Use the volecity verlet to update"<<std::endl;
                    q += params_.timeStep*v;
                    computeForceAndHessian(q, oldq, F, H);
                    v += params_.timeStep*Minv_*F;
                }
                break;
            }
            default:
                break;
        }
        Eigen::VectorXd curq = q;
        int nRods = getNumRigidRods();
        Eigen::VectorXd x(q.size() + nRods);
        x.setZero();
        x.topRows(q.size()) = q;
        
        Eigen::VectorXd x0 = x;
        
        flag = newtonSolver(x, [this, x0](Eigen::VectorXd x, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF)
                     {
                         this->computeStepProjection(x, x0, F, gradF);
                     });
        if(flag)    // successfully update
        {
            q = x.topRows(q.size());
            v += (q-curq)/params_.timeStep;
        }
        else
        {
            q = oldq;
            v = oldv;
            std::cout<<"Failed to update, please use a smaller time step, current time step is: "<<params_.timeStep<<std::endl;
        }
    }
    else if(params_.constraintHandling == params_.CH_LAGRANGEMULT)
    {
        q += params_.timeStep*v;
        computeForceAndHessian(q, oldq, F, H);
        Eigen::SparseMatrix<double> gradG;
        Eigen::VectorXd g,lambda;
        computeContraintsAndGradient(q, g, &gradG);
        
        int nRods = getNumRigidRods();
        lambda.resize(nRods);
        lambda.setZero();
        
       
        flag = newtonSolver(lambda, [this, q, v, F](Eigen::VectorXd lambda, Eigen::VectorXd &force, Eigen::SparseMatrix<double> *gradF)
                            {
                                this->computeLagrangeMultiple(lambda, q, v, F, force, gradF);
                            });
        if(flag)
        {
            Eigen::VectorXd g;
            Eigen::SparseMatrix<double> gradG;
            computeContraintsAndGradient(q, g, &gradG);
            v += params_.timeStep*Minv_*(F + gradG.transpose()*lambda);
        }
        else
        {
            q = oldq;
            v = oldv;
            std::cout<<"Failed to update, please use a smaller time step, current time step is: "<<params_.timeStep<<std::endl;
        }
    }
    return flag;
}

void GooHook::updatebyVelocityVerletUnconstraint(Eigen::VectorXd & q, Eigen::VectorXd & v)
{
	SparseMatrix<double> H;
	VectorXd oldq = q;
	VectorXd F; 
	q += params_.timeStep*v;
	computeForceAndHessian(q, oldq, F, H);
    std::vector<Eigen::Triplet<double>> penaltyCoeffs;
    processPenaltyForce(q, oldq, F, penaltyCoeffs);
	v += params_.timeStep*Minv_*F;

}

void GooHook::updatebyImpliciyMidpointUnconstraint(Eigen::VectorXd & q, Eigen::VectorXd & v, const Eigen::VectorXd prevq)
{
	Eigen::VectorXd currentq = q;
	Eigen::SparseMatrix<double> idMat(q.size(), q.size());
	idMat.setIdentity();
	

	for (int i = 0; i < params_.NewtonMaxIters; i++)
	{
		Eigen::VectorXd force;
		Eigen::SparseMatrix<double> H;
		
		// H = -gradF
		computeForceAndHessian((currentq + q) / 2, (q + prevq) / 2, force, H);
        std::vector<Eigen::Triplet<double>> penaltyCoeffs;
        processPenaltyForce((currentq + q) / 2, (q + prevq) / 2, force, penaltyCoeffs);
        SparseMatrix<double> penaltyMat(q.size(), q.size());
        penaltyMat.setFromTriplets(penaltyCoeffs.begin(), penaltyCoeffs.end());
        H += penaltyMat;

		Eigen::VectorXd fVal = currentq - q - params_.timeStep*v - 0.5*params_.timeStep*params_.timeStep*Minv_*force;


		//std::cout << "Feval: " << fVal.norm() << std::endl;
		if (fVal.norm() < params_.NewtonTolerance)
			break;

		Eigen::SparseMatrix<double> gradFVal = idMat - params_.timeStep * params_.timeStep * Minv_ * (-H) / 4.0;

		//Warning: gradF is unsymmetry anymore.
		Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver; 
		solver.compute(gradFVal);

		Eigen::VectorXd deltaq = solver.solve(-fVal);
		//std::cout <<"Deltaq: "<< deltaq.norm() << std::endl;
		currentq = currentq + deltaq;
	}

	v = 2.0 * (currentq - q) / params_.timeStep - v;
	q = currentq;
}

void GooHook::computeLagrangeMultiple(Eigen::VectorXd lambda, Eigen::VectorXd pos, Eigen::VectorXd vel, Eigen::VectorXd regularForce, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF)
{
    Eigen::VectorXd newq, newv;
    Eigen::VectorXd g;
    Eigen::SparseMatrix<double> gradG;
    computeContraintsAndGradient(pos, g, &gradG);
    newv = vel + params_.timeStep*Minv_*(regularForce + gradG.transpose()*lambda);
    newq = pos + params_.timeStep*newv;
    if(gradF)
    {
        computeContraintsAndGradient(newq, F, gradF);
        Eigen::SparseMatrix<double> W;
        W = params_.timeStep * params_.timeStep * Minv_ *gradG.transpose();
        *gradF = (*gradF) * W;
    }
}

void GooHook::computeContraintsAndGradient(const Eigen::VectorXd q, Eigen::VectorXd &g, Eigen::SparseMatrix<double> *gradG)
{
    int nConnectors = connectors_.size();
    int nRods = getNumRigidRods();
    g.resize(nRods);
    g.setZero();
    
    std::vector<Eigen::Triplet<double> > dgTriplet;
    int rodIdx = 0;
    
    for(int i=0;i<nConnectors;i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
            continue;
        auto rod = *(RigidRod *)connectors_[i];
        Eigen::Vector2d p1 = q.segment<2>(2*rod.p1);
        Eigen::Vector2d p2 = q.segment<2>(2*rod.p2);
        
        g(rodIdx) = (p1-p2).squaredNorm() - rod.length*rod.length;
        
        if(gradG)
        {
            Eigen::Vector2d localF = 2*(p1-p2);
            
            dgTriplet.push_back(Eigen::Triplet<double>(rodIdx, 2*rod.p1, localF(0)));
            dgTriplet.push_back(Eigen::Triplet<double>(rodIdx, 2*rod.p1 + 1, localF(1)));
            dgTriplet.push_back(Eigen::Triplet<double>(rodIdx, 2*rod.p2, -localF(0)));
            dgTriplet.push_back(Eigen::Triplet<double>(rodIdx, 2*rod.p2 + 1, -localF(1)));
        }
        rodIdx ++;
    }
    if(gradG)
    {
        gradG->resize(nRods, q.size());
        gradG->setFromTriplets(dgTriplet.begin(), dgTriplet.end());
    }
}


void GooHook::computeStepProjection(Eigen::VectorXd x, Eigen::VectorXd x0, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF)
{
    int rodIdx = 0;
    int nParticles = particles_.size();
    int nConnectors = connectors_.size();
    int nRods = x0.size() - 2*nParticles;
    Eigen::VectorXd g(nRods);
    Eigen::VectorXd f(2*nParticles), dg(2*nParticles);
    f = x.topRows(2*nParticles) - x0.topRows(2*nParticles);
    dg.setZero();
    
    std::vector<Eigen::Triplet<double> > mInvTriplet, liftTriplet, IdgTriplet;
    
    for(int i=0; i<nParticles; i++)
    {
        liftTriplet.push_back(Eigen::Triplet<double>(2*i, 2*i, 1));
        liftTriplet.push_back(Eigen::Triplet<double>(2*i + 1, 2*i + 1, 1));
       
        IdgTriplet.push_back(Eigen::Triplet<double>(2*i, 2*i, 1));
        IdgTriplet.push_back(Eigen::Triplet<double>(2*i + 1, 2*i + 1, 1));
    }
    
    for(int i = 0; i<nConnectors; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
            continue;
        auto rod = *(RigidRod *)connectors_[i];
        Eigen::Vector2d p1 = x.segment<2>(2*rod.p1);
        Eigen::Vector2d p2 = x.segment<2>(2*rod.p2);
        
        Eigen::Vector2d localF = 2*(p1-p2);
        
        mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p1, 2*nParticles + rodIdx, localF(0)));
        mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p1 + 1, 2*nParticles + rodIdx, localF(1)));
        mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p2, 2*nParticles + rodIdx, -localF(0)));
        mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p2 + 1, 2*nParticles + rodIdx, -localF(1)));
        
        IdgTriplet.push_back(Eigen::Triplet<double>(2*nParticles + rodIdx, 2*rod.p1, localF(0)));
        IdgTriplet.push_back(Eigen::Triplet<double>(2*nParticles + rodIdx, 2*rod.p1 + 1, localF(1)));
        IdgTriplet.push_back(Eigen::Triplet<double>(2*nParticles + rodIdx, 2*rod.p2, -localF(0)));
        IdgTriplet.push_back(Eigen::Triplet<double>(2*nParticles + rodIdx, 2*rod.p2 + 1, -localF(1)));
        
        
        dg.segment<2>(2*rod.p1) += x(2*nParticles + rodIdx) * localF;
        dg.segment<2>(2*rod.p2) += -x(2*nParticles + rodIdx) * localF;
        
        g(rodIdx) = (p1-p2).squaredNorm() - rod.length*rod.length;
        
        Eigen::Matrix2d localH;
        localH.setIdentity();
        localH = 2*localH;
        
        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p1+j, 2*rod.p1+k, x(2*nParticles + rodIdx)* localH.coeff(j,k)));
                mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p2+j, 2*rod.p2+k, x(2*nParticles + rodIdx) * localH.coeff(j,k)));
                mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p1+j, 2*rod.p2+k, -x(2*nParticles + rodIdx)*localH.coeff(j,k)));
                mInvTriplet.push_back(Eigen::Triplet<double>(2*rod.p2+j, 2*rod.p1+k, -x(2*nParticles + rodIdx)*localH.coeff(j,k)));
            }
        
        
        rodIdx++;
    }
//    std::cout<<"lambda*dg"<<std::endl;
//    std::cout<<dg<<std::endl<<std::endl;
//
//    std::cout<<"f"<<std::endl;
    f += Minv_ * dg;
//    std::cout<<f<<std::endl;
    
    F.resize(x.size());
    F.topRows(2*nParticles) = f;
    F.bottomRows(nRods) = g;
    
    if(gradF == NULL)
        return;
    
    gradF->resize(x.size(), x.size());
    gradF->setZero();
    
    SparseMatrix<double> liftMat(x.size(), 2*nParticles);
    liftMat.setFromTriplets(liftTriplet.begin(), liftTriplet.end());
    
    SparseMatrix<double> IdgMat(x.size(),x.size());
    IdgMat.setFromTriplets(IdgTriplet.begin(), IdgTriplet.end());
    
    SparseMatrix<double> mInvMat(2*nParticles, x.size());
    mInvMat.setFromTriplets(mInvTriplet.begin(), mInvTriplet.end());

//    std::cout<<"[lambda * H, dg]"<<std::endl;
//    std::cout<<mInvMat.toDense()<<std::endl<<std::endl;
    *gradF = IdgMat + liftMat * Minv_ * mInvMat;
//    std::cout<<"Grad F"<<std::endl;
//    std::cout<<gradF.toDense()<<std::endl<<std::endl;
    
    
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

void GooHook::computeMassInverse(Eigen::SparseMatrix<double> &Minv)
{
    int ndofs = 2*int(particles_.size());

    Minv.resize(ndofs, ndofs);
    Minv.setZero();

    std::vector<Eigen::Triplet<double> > Minvcoeffs;
    for(int i=0; i<ndofs/2; i++)
    {
        Minvcoeffs.push_back(Eigen::Triplet<double>(2*i,   2*i,   1.0/getTotalParticleMass(i)));
        Minvcoeffs.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, 1.0/getTotalParticleMass(i)));
    }

    Minv.setFromTriplets(Minvcoeffs.begin(), Minvcoeffs.end());
}

void GooHook::computeForceAndHessian(const VectorXd &q, const VectorXd &qprev, Eigen::VectorXd &F, SparseMatrix<double> &H)
{
    F.resize(q.size());
    F.setZero();
    H.resize(q.size(), q.size());
    H.setZero();

    std::vector<Eigen::Triplet<double> > Hcoeffs;
    if(params_.gravityEnabled)
        processGravityForce(F);
    if(params_.springsEnabled)
        processSpringForce(q, F, Hcoeffs);
    if(params_.dampingEnabled)
        processDampingForce(q, qprev, F, Hcoeffs);
    if(params_.floorEnabled)
        processFloorForce(q, qprev, F, Hcoeffs);
    if(params_.bendingEnabled)
        processBendingForce(q, qprev, F, Hcoeffs);
    H.setFromTriplets(Hcoeffs.begin(), Hcoeffs.end());
}

void GooHook::processGravityForce(VectorXd &F)
{
    int nparticles = (int)particles_.size();
    for(int i=0; i<nparticles; i++)
    {
        if(!particles_[i].fixed)
        {
            F[2*i+1] += params_.gravityG*getTotalParticleMass(i);
        }
    }
}

void GooHook::processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nsprings = (int)connectors_.size();

    for(int i=0; i<nsprings; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring &s = *(Spring *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*s.p1);
        Vector2d p2 = q.segment<2>(2*s.p2);
        double dist = (p2-p1).norm();
        Vector2d localF = s.stiffness*(dist-s.restlen)/dist * (p2-p1);
        F.segment<2>(2*s.p1) += localF;
        F.segment<2>(2*s.p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = s.stiffness * (1.0 - s.restlen/dist)*I;
        localH += s.stiffness*s.restlen*(p2-p1)*(p2-p1).transpose()/dist/dist/dist;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p1+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p2+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p2+k, -localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p1+k, -localH.coeff(j,k)));
            }
    }
}

void GooHook::processDampingForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nsprings = (int)connectors_.size();

    for(int i=0; i<nsprings; i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;
        Spring &s = *(Spring *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*s.p1);
        Vector2d p2 = q.segment<2>(2*s.p2);
        Vector2d p1prev = qprev.segment<2>(2*s.p1);
        Vector2d p2prev = qprev.segment<2>(2*s.p2);

        Vector2d relvel = (p2 - p2prev)/params_.timeStep - (p1 - p1prev)/params_.timeStep;
        Vector2d localF = params_.dampingStiffness*relvel;
        F.segment<2>(2*s.p1) += localF;
        F.segment<2>(2*s.p2) -= localF;

        Matrix2d I;
        I << 1, 0, 0, 1;
        Matrix2d localH = params_.dampingStiffness*I/params_.timeStep;

        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p1+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p2+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p2+k, -localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p1+k, -localH.coeff(j,k)));
            }
    }
}

void GooHook::processFloorForce(const VectorXd &q, const VectorXd &qprev, VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nparticles = particles_.size();

    double basestiffness = 10000;
    double basedrag = 1000.0;
	double baseradius = 0.02;

    for(int i=0; i<nparticles; i++)
    {
		// Make the simulation more makes sense
		double radius = baseradius*sqrt(getTotalParticleMass(i));
        if(q[2*i+1] < -0.5+radius && ! particles_[i].fixed)
        {
            double vel = (q[2*i+1]-qprev[2*i+1])/params_.timeStep;
            double dist = -0.5+ radius - q[2*i+1];

            F[2*i+1] += basestiffness*dist - basedrag*dist*vel;

            H.push_back(Eigen::Triplet<double>(2*i+1, 2*i+1, basestiffness
                - 0.5*basedrag/params_.timeStep
                + basedrag*qprev[2*i+1]/params_.timeStep
                - 2.0*basedrag*q[2*i+1]/params_.timeStep));
        }
    }
}

void GooHook::processPenaltyForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H)
{
    int nConnects = connectors_.size();
    
    for(int i=0;i<nConnects;i++)
    {
        if(connectors_[i]->getType() != SimParameters::CT_RIGIDROD)
            continue;
        RigidRod &s = *(RigidRod *)connectors_[i];
        Vector2d p1 = q.segment<2>(2*s.p1);
        Vector2d p2 = q.segment<2>(2*s.p2);
        double squaredDist = (p2-p1).squaredNorm();
        Vector2d localF = -4 * params_.penaltyStiffness * (squaredDist - s.length*s.length)*(p1-p2);
        F.segment<2>(2*s.p1) += localF;
        F.segment<2>(2*s.p2) -= localF;
        
        Matrix2d I;
        I << 1,0,0,1;
        
        // H = -dF
        Eigen::Matrix2d localH = 4*params_.penaltyStiffness*(squaredDist-s.length*s.length)*I + 8*params_.penaltyStiffness*(p1-p2)*(p1-p2).transpose();
        
        for(int j=0; j<2; j++)
            for(int k=0; k<2;k++)
            {
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p1+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p2+k, localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p1+j, 2*s.p2+k, -localH.coeff(j,k)));
                H.push_back(Eigen::Triplet<double>(2*s.p2+j, 2*s.p1+k, -localH.coeff(j,k)));
            }
        
        
    }
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

double GooHook::ptSegmentDist(const Vector2d &p, const Vector2d &q1, const Vector2d &q2)
{
    double t = (p-q1).dot(q2-q1) / (q2-q1).dot(q2-q1);
    double linedistsq = (q1 + t*(q2-q1) - p).squaredNorm();
    double q1dist = (p-q1).squaredNorm();
    double q2dist = (p-q2).squaredNorm();
    double mindistsq = std::min(linedistsq, std::min(q1dist, q2dist));
    return sqrt(mindistsq);
}

void GooHook::detectSawedConnectors(std::set<int> &connectorsToDelete)
{
    for(int i=0; i<(int)connectors_.size(); i++)
    {
        Vector2d pos1 = particles_[connectors_[i]->p1].pos;
        Vector2d pos2 = particles_[connectors_[i]->p2].pos;
        double maxx = std::max(pos1[0], pos2[0]);
        double minx = std::min(pos1[0], pos2[0]);
        double maxy = std::max(pos1[1], pos2[1]);
        double miny = std::min(pos1[1], pos2[1]);
        for(std::vector<Saw>::iterator saw = saws_.begin(); saw != saws_.end(); ++saw)
        {
            Vector2d sawpos = saw->pos;
            double sawr = saw->radius;

            if(sawpos[0] - sawr > maxx || sawpos[0] + sawr < minx || sawpos[1] - sawr > maxy || sawpos[1] + sawr < miny)
                continue;

            double sawspringdist = ptSegmentDist(sawpos, pos1, pos2);
            if(sawspringdist <= sawr)
            {
                connectorsToDelete.insert(i);
                break;
            }
        }
    }
}

void GooHook::detectSawedParticles(std::set<int> &particlesToDelete)
{
    for(int i=0; i<(int)particles_.size(); i++)
    {
        Vector2d partpos = particles_[i].pos;

        if(!(fabs(partpos[0]) < 2 && fabs(partpos[1]) < 2))
        {
            particlesToDelete.insert(i);
            break;
        }

        for(std::vector<Saw>::iterator it = saws_.begin(); it != saws_.end(); ++it)
        {
            Vector2d sawpos = it->pos;
            double sqdist = (sawpos-partpos).squaredNorm();
            if(sqdist < it->radius*it->radius)
            {
                particlesToDelete.insert(i);
                break;
            }
        }
    }
}

void GooHook::deleteSawedObjects()
{
    std::set<int> particlestodelete;
    std::set<int> connectorstodelete;
    std::set<int> bendingtodelete;
    detectSawedParticles(particlestodelete);
    detectSawedConnectors(connectorstodelete);

    std::vector<Particle, Eigen::aligned_allocator<Particle>> newparticles;
    std::vector<Connector *> newconnectors;
    std::vector<BendingStencil> newbending;
    std::vector<int> remainingparticlemap;
    std::vector<int> remainingbendingmap;

    if(!particlestodelete.empty())
    {
        for(int i=0; i<(int)connectors_.size(); i++)
        {
            if(particlestodelete.count(connectors_[i]->p1) || particlestodelete.count(connectors_[i]->p2))
                connectorstodelete.insert(i);
        }

        for(int i=0; i<(int)particles_.size(); i++)
        {
            if(particlestodelete.count(i) == 0)
            {
                remainingparticlemap.push_back(newparticles.size());
                newparticles.push_back(particles_[i]);
            }
            else
                remainingparticlemap.push_back(-1);
        }
    }
    if(!connectorstodelete.empty())
    {
        for(std::set<int>::iterator it = connectorstodelete.begin(); it != connectorstodelete.end(); ++it)
        {
            for(std::set<int>::iterator bit = connectors_[*it]->associatedBendingStencils.begin(); bit != connectors_[*it]->associatedBendingStencils.end(); ++bit)
            {
                bendingtodelete.insert(*bit);
            }
        }
        for(int i=0; i<(int)connectors_.size(); i++)
        {
            if(connectorstodelete.count(i) == 0)
            {
                newconnectors.push_back(connectors_[i]);
            }
            else
                delete connectors_[i];
        }
    }
    if(!bendingtodelete.empty())
    {
        int newidx=0;
        for(int i=0; i<(int)bendingStencils_.size(); i++)
        {
            if(bendingtodelete.count(i) == 0)
            {
                newbending.push_back(bendingStencils_[i]);
                remainingbendingmap.push_back(newidx++);
            }
            else
                remainingbendingmap.push_back(-1);
        }
    }

    if (!connectorstodelete.empty() || !particlestodelete.empty())
    {
        if (!connectorstodelete.empty())
            connectors_ = newconnectors;
        if (!bendingtodelete.empty())
        {
            bendingStencils_ = newbending;
            for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
            {
                std::set<int> newass;
                for (std::set<int>::iterator sit = (*it)->associatedBendingStencils.begin(); sit != (*it)->associatedBendingStencils.end(); ++sit)
                {
                    if (bendingtodelete.count(*sit) == 0)
                        newass.insert(remainingbendingmap[*sit]);
                }
                (*it)->associatedBendingStencils = newass;
            }
        }
        if (!particlestodelete.empty())
        {
            particles_ = newparticles;
            for (std::vector<Connector *>::iterator it = connectors_.begin(); it != connectors_.end(); ++it)
            {
                (*it)->p1 = remainingparticlemap[(*it)->p1];
                (*it)->p2 = remainingparticlemap[(*it)->p2];
            }
            for (std::vector<BendingStencil>::iterator it = bendingStencils_.begin(); it != bendingStencils_.end(); ++it)
            {
                it->p1 = remainingparticlemap[it->p1];
                it->p2 = remainingparticlemap[it->p2];
                it->p3 = remainingparticlemap[it->p3];
            }
        }
    }
}

void GooHook::pruneOverstrainedSprings()
{
    int nsprings = connectors_.size();

    std::vector<int> toremove;
    for (int i = 0; i < nsprings; i++)
    {
        if (connectors_[i]->getType() != SimParameters::CT_SPRING)
            continue;

        Spring &s = *(Spring *)connectors_[i];
        if (s.canSnap)
        {
            Vector2d srcpos = particles_[s.p1].pos;
            Vector2d dstpos = particles_[s.p2].pos;
            double dist = (dstpos - srcpos).norm();

            double strain = (dist - s.restlen) / s.restlen;
            if (strain > params_.maxSpringStrain)
                toremove.push_back(i);
        }
    }

    for (std::vector<int>::reverse_iterator it = toremove.rbegin(); it != toremove.rend(); ++it)
    {
        assert(connectors_[*it]->associatedBendingStencils.empty());
        delete connectors_[*it];
        connectors_.erase(connectors_.begin() + *it);
    }
}




//////////////////////////////////////////////////////////////////////////////////////
////    Helper Function
/////////////////////////////////////////////////////////////////////////////////////

void GooHook::testForceDifferential()
{

	Eigen::VectorXd q, vel, qPrev;
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
	}
}

void GooHook::testStepProjection()
{
    Eigen::VectorXd q, vel, qPrev;
    buildConfiguration(q, vel, qPrev);
    qPrev.resize(q.size());
    for (int i = 0; i<(int)particles_.size(); i++)
    {
        qPrev.segment<2>(2 * i) = particles_[i].prevpos;
    }
    

    int nRods = getNumRigidRods();


    Eigen::VectorXd x(q.size() + nRods);
    x.setOnes();
    x.topRows(q.size()) = q;
    
    Eigen::VectorXd x0 = x;
    
    Eigen::SparseMatrix<double> gradF(x.size(), x.size());
    gradF.setZero();
    Eigen::VectorXd f = Eigen::VectorXd::Zero(x.size());

    computeStepProjection(x, x0, f, &gradF);

    Eigen::VectorXd direction = Eigen::VectorXd::Random(x.size());
    
    for(int i=0;i<particles_.size();i++)
    {
        if(particles_[i].fixed)
            direction.segment(2*i, 2).setZero();
    }
    

    direction.normalize();
    std::cout<<direction<<std::endl;
    for (int k = 1; k <= 12; k++)
    {
        double eps = pow(10, -k);
        
        VectorXd epsF = Eigen::VectorXd::Zero(x.size());
        Eigen::SparseMatrix<double>  epsgradF(x.size(), x.size());
        computeStepProjection(x+eps*direction, x0, epsF, &epsgradF);
        
        std::cout << "EPS is: " << eps << std::endl;
        std::cout << "Norm of Finite Difference is: " << (epsF - f).norm() / eps << std::endl;
        std::cout << "Norm of Directinal Gradient is: " << (gradF*direction).norm() << std::endl;
        std::cout << "The difference between above two is: " << ((epsF - f) / eps - gradF*direction).norm() << std::endl << std::endl;

    }
}

void GooHook::testLagrangeMultiple()
{
    Eigen::VectorXd q, vel, qPrev;
    buildConfiguration(q, vel, qPrev);
    qPrev.resize(q.size());
    for (int i = 0; i<(int)particles_.size(); i++)
    {
        qPrev.segment<2>(2 * i) = particles_[i].prevpos;
    }
    
    Eigen::VectorXd F;
    Eigen::SparseMatrix<double> H;
    
    q += params_.timeStep * vel;
    computeForceAndHessian(q, qPrev, F, H);
    int nRods = getNumRigidRods();
    
    
    Eigen::VectorXd lambda = Eigen::VectorXd::Random(nRods);
    
    Eigen::SparseMatrix<double> gradF(lambda.size(), lambda.size());
    gradF.setZero();
    Eigen::VectorXd f = Eigen::VectorXd::Zero(lambda.size());
    
    computeLagrangeMultiple(lambda, q, vel, F, f, &gradF);
    
    Eigen::VectorXd direction = Eigen::VectorXd::Random(lambda.size());
    
    
    direction.normalize();
    std::cout<<direction<<std::endl;
    for (int k = 1; k <= 12; k++)
    {
        double eps = pow(10, -k);
        
        VectorXd epsF = Eigen::VectorXd::Zero(lambda.size());
        Eigen::SparseMatrix<double>  epsgradF(lambda.size(), lambda.size());
        computeLagrangeMultiple(lambda+eps*direction, q, vel, F, epsF, &epsgradF);
        
        std::cout << "EPS is: " << eps << std::endl;
        std::cout << "Norm of Finite Difference is: " << (epsF - f).norm() / eps << std::endl;
        std::cout << "Norm of Directinal Gradient is: " << (gradF*direction).norm() << std::endl;
        std::cout << "The difference between above two is: " << ((epsF - f) / eps - gradF*direction).norm() << std::endl << std::endl;
        
    }
}

void GooHook::testConstriantAndGradient()
{
    Eigen::VectorXd q, vel, qPrev;
    buildConfiguration(q, vel, qPrev);
    qPrev.resize(q.size());
    for (int i = 0; i<(int)particles_.size(); i++)
    {
        qPrev.segment<2>(2 * i) = particles_[i].prevpos;
    }
    
    
    int nRods = getNumRigidRods();
    
    Eigen::VectorXd g;
    Eigen::SparseMatrix<double> gradG;
    computeContraintsAndGradient(q, g, &gradG);
    
    Eigen::VectorXd direction = Eigen::VectorXd::Random(q.size());
    
    for(int i=0;i<particles_.size();i++)
    {
        if(particles_[i].fixed)
            direction.segment(2*i, 2).setZero();
    }
    
    
    direction.normalize();
    std::cout<<direction<<std::endl;
    for (int k = 1; k <= 12; k++)
    {
        double eps = pow(10, -k);
        
        VectorXd epsF;
        Eigen::SparseMatrix<double>  epsgradF;
        computeContraintsAndGradient(q+eps*direction, epsF, &epsgradF);
        
        std::cout << "EPS is: " << eps << std::endl;
        std::cout << "Norm of Finite Difference is: " << (epsF - g).norm() / eps << std::endl;
        std::cout << "Norm of Directinal Gradient is: " << (gradG*direction).norm() << std::endl;
        std::cout << "The difference between above two is: " << ((epsF - g) / eps - gradG*direction).norm() << std::endl << std::endl;
        
    }
}

void GooHook::saveConfiguration(std::string filePath)
{
	int nParticles = particles_.size();
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
	outfile.close();
}

void GooHook::loadConfiguration(std::string filePath)
{
	std::ifstream infile(filePath);
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
	}
}
