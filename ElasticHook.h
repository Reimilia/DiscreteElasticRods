#ifndef ELASTICHOOK_H
#define ELASTICHOOK_H

#include <deque>

#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include "PhysicsHook.h"
#include "SceneObjects.h"
#include "SimParameters.h"
#include "VectorMath.h"
#include "ElasticRod.h"
#include "RigidBodyInstance.h"
#include "RigidBodyTemplate.h"

struct MouseClick
{
	double x;
	double y;
	SimParameters::ClickMode mode;
};

class ElasticHook : public PhysicsHook
{
public:
	ElasticHook() : PhysicsHook() { launch_ = false; ballTemplate_ = NULL; rodTemplate_ = NULL;  isNewRod_ = true; };

	virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

	virtual void initSimulation();

	virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir, int button);
	//virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer, int button);
	//virtual bool mouseMoved(igl::opengl::glfw::Viewer &viewer, int button);

	virtual void updateRenderGeometry();

	virtual void tick();

	virtual bool simulateOneStep();

	virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
	{
		viewer.data().clear();
		viewer.data().set_mesh(renderQ, renderF);
	}

	void saveConfiguration(std::string filePath);
	void loadConfiguration(std::string filePath);


public:



private:
	bool launch_;
	Eigen::Vector3d launchPos_;
	Eigen::Vector3d launchDir_;

	bool isNewRod_;

	// Load sphere and rod model
	// These two elements are only for rendering
	// So we will not assign them any rigid body property
	RigidBodyTemplate *ballTemplate_;
	RigidBodyTemplate *rodTemplate_;

	void loadMesh();

	SimParameters params_;
	double time_;

	//std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;

	//Rigid Body Elements
    std::vector<RigidBodyTemplate *> templates_;
	std::vector<RigidBodyInstance *> bodies_;

	//Elastic Rods Elements
	//Note that the model assumes only connecting points are d.o.f.
	//So I use a lazy way to render the rods_
	std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;
	std::vector<ElasticRod *> rods_;

	std::mutex message_mutex;
	std::deque<MouseClick> mouseClicks_;

	Eigen::MatrixXd renderQ;
	Eigen::MatrixXi renderF;
	Eigen::MatrixXd renderC;

	Eigen::SparseMatrix<double> mInv_;

	void addParticle(double x, double y, double z);

	void buildConfiguration(Eigen::VectorXd &, Eigen::VectorXd &);
	void unbuildConfiguration(const Eigen::VectorXd &, const Eigen::VectorXd &);

	void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
	bool numericalIntegration();

	//void updatebyVelocityVerletUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v);
	//void updatebyImpliciyMidpointUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v, const Eigen::VectorXd prevq);

	void computeForce(Eigen::VectorXd &F);
	//void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);

	void computeContraintsAndGradient(const Eigen::VectorXd &, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &gradF);
	void computeLagrangeMultiple(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::SparseMatrix<double> &, const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &);

	void computeGravityForce(Eigen::VectorXd &F);
	//void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
	
	bool newtonSolver(Eigen::VectorXd &x, std::function<void(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> _computeForceAndHessian);

	void testProcess();
};

#endif