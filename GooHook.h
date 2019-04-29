#include "PhysicsHook.h"
#include "SceneObjects.h"
#include <deque>
#include "SimParameters.h"
#include "VectorMath.h"
#include <Eigen/Sparse>
#include <Eigen/StdVector>

struct MouseClick
{
    double x;
    double y;
    SimParameters::ClickMode mode;
};

class GooHook : public PhysicsHook
{
public:
	GooHook() : PhysicsHook() {};

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation();

    virtual void mouseClicked(double x, double y, int button)
    {
        message_mutex.lock();
        {
            MouseClick mc;
            mc.x = x;
            mc.y = y;
            mc.mode = params_.clickMode;
            mouseClicks_.push_back(mc);
        }
        message_mutex.unlock();
    }

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
        //viewer.data().set_colors(renderC);
    }

	// Added Test function
	void testForceDifferential();

	void saveConfiguration(std::string filePath);
	void loadConfiguration(std::string filePath);
    
    
public:
    


private:
	// Load sphere and rod model
	Eigen::MatrixXd ballV;
	Eigen::MatrixXi ballF;
	Eigen::MatrixXd rodV;
	Eigen::MatrixXi rodF;

	void loadMesh();

    SimParameters params_;
    double time_;

	//TODO: replace each conncector_ as rigid body template
	//Note each connector for elastic rod will have its own bishop frame.
    std::vector<Particle, Eigen::aligned_allocator<Particle> > particles_;
    std::vector<Connector *> connectors_;
    std::vector<BendingStencil> bendingStencils_;


    std::mutex message_mutex;
    std::deque<MouseClick> mouseClicks_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd renderC;
    
    //Eigen::SparseMatrix<double> Minv_;

    void addParticle(double x, double y);
    double getTotalParticleMass(int idx);
    int getNumRigidRods();

	void updateQuasiStaticFrame();
	void updateBishopFrame();


	void buildConfiguration(Eigen::VectorXd &q, Eigen::VectorXd &v, Eigen::VectorXd &qprev);
    void unbuildConfiguration(const Eigen::VectorXd &q, const Eigen::VectorXd &v);

    //void computeMassInverse(Eigen::SparseMatrix<double> &Minv);
	//void computeMass(Eigen::SparseMatrix<double> &M);
    bool numericalIntegration(Eigen::VectorXd &q, Eigen::VectorXd &v, Eigen::VectorXd &prevq);

	//void updatebyVelocityVerletUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v);
	//void updatebyImpliciyMidpointUnconstraint(Eigen::VectorXd &q, Eigen::VectorXd &v, const Eigen::VectorXd prevq);

    void computeForceAndHessian(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, Eigen::SparseMatrix<double> &H);
    
    //void computeContraintsAndGradient(const Eigen::VectorXd q, Eigen::VectorXd &F, Eigen::SparseMatrix<double> *gradF);
    
    void processGravityForce(Eigen::VectorXd &F);
    void processSpringForce(const Eigen::VectorXd &q, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
    void processDampingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
    void processFloorForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
    void processPenaltyForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
	void processBendingForce(const Eigen::VectorXd &q, const Eigen::VectorXd &qprev, Eigen::VectorXd &F, std::vector<Eigen::Triplet<double> > &H);
    
    bool newtonSolver(Eigen::VectorXd &x, std::function<void(Eigen::VectorXd, Eigen::VectorXd &, Eigen::SparseMatrix<double> *)> _computeForceAndHessian);
    bool takeOneStep(double &ratio, Eigen::VectorXd &x, std::function<void(Eigen::VectorXd, Eigen::VectorXd &, Eigen::SparseMatrix<double> *)> _computeForceAndHessian);

};
