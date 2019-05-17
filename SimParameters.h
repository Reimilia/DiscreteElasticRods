#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

#include <Eigen/Core>

struct SimParameters
{
    SimParameters()
    {
        constraintHandling = CH_PENALTY;
		integrator = TI_VELOCITY_VERLET;
        timeStep = 0.001;
        NewtonMaxIters = 20;
        NewtonTolerance = 1e-10;
        penaltyStiffness = 1e5;

        gravityEnabled = true;
        gravityG = -9.8;
        //springsEnabled = true;
        //springStiffness = 100;
        //maxSpringStrain = 0.2;
        //dampingEnabled = true;
        //dampingStiffness = 1.0;
        floorEnabled = true;
        //bendingEnabled = true;
		//twistEnabled = true;

        clickMode = CM_ADDPARTICLE;
        connectorType = CT_FLEXROD;
        particleMass = 0.0;
        //maxSpringDist = 1e5;
        particleFixed = false;

        rodDensity = 1;
		rodTwistStiffness = 0.075;
		//rodStretchingStiffness = 100;
        //rodBendingStiffness = 0.05;
        rodSegments = 5;

        //sawRadius= 0.1;
		rodBendingModulus << 0.1, -0.03, -0.03, 0.05;
		boundaryCondition = BC_FREE;

		leftAngularVelocity = 0;
		rightAngularVelocity = M_PI;
    }

	enum BoundaryCondition
	{
		BC_FREE,
		BC_FIXED_END,
		BC_RIGIDBODY_END
	};
    enum ClickMode {CM_ADDPARTICLE, CM_FREE};

    enum ConstraintHandling {CH_PENALTY, CH_STEPPROJECT, CH_LAGRANGEMULT};
	enum TimeIntegrator {TI_VELOCITY_VERLET, TI_IMPLICIT_MIDPOINT};
    enum ConnectorType {CT_FLEXROD, CT_SPRING, CT_RIGIDROD, CT_ELASTICROD};

	BoundaryCondition boundaryCondition;
    ConstraintHandling constraintHandling;
	TimeIntegrator integrator;
    double timeStep;
    double NewtonTolerance;
    int NewtonMaxIters;
    double penaltyStiffness;

    bool gravityEnabled;
    double gravityG;
    //bool springsEnabled;
    //bool bendingEnabled;
	//bool twistEnabled;
    //double springStiffness;
    //double maxSpringStrain;
    bool floorEnabled;
    //bool dampingEnabled;
    //double dampingStiffness;

    ClickMode clickMode;
    ConnectorType connectorType;
    double particleMass;
    bool particleFixed;
    //double sawRadius;

	double rodDensity;
	double rodTwistStiffness;
	int rodSegments;
	Eigen::Matrix2d rodBendingModulus;
	double leftAngularVelocity, rightAngularVelocity;
};

#endif