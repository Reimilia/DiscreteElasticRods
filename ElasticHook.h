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
	ElasticHook() : PhysicsHook() {};

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

	// Added Test function
	void testForceDifferential();

	void saveConfiguration(std::string filePath);
	void loadConfiguration(std::string filePath);


public:



private:
	bool launch_;
	Eigen::Vector3d launchPos_;
	Eigen::Vector3d launchDir_;


	// Load sphere and rod model
	// These two elements are only for rendering
	// So we will not assign them any rigid body property
	Eigen::MatrixXd ballV;
	Eigen::MatrixXi ballF;
	Eigen::MatrixXd rodV;
	Eigen::MatrixXi rodF;

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
	std::vector<ElasticRod *> rods_;

	std::mutex message_mutex;
	std::deque<MouseClick> mouseClicks_;

	Eigen::MatrixXd renderQ;
	Eigen::MatrixXi renderF;
	Eigen::MatrixXd renderC;

	//Eigen::SparseMatrix<double> Minv_;

	void addParticle(double x, double y, double z);

	//void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
	//void computeMass(Eigen::SparseMatrix<double> &M);
	bool numericalIntegration();

	//void updatebyVelocityVerletUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v);
	//void updatebyImpliciyMidpointUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v, const Eigen::VectorXd prevq);

	void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);

	//void computeContraintsAndGradient(const Eigen::VectorXd q, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF);

	void processGravityForce(Eigen::VectorXd &F);
	void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
	
	bool newtonSolver(Eigen::VectorXd &x, std::function<void(Eigen::VectorXd &, Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> _computeForceAndHessian);

};

#endif