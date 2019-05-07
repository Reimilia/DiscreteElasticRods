#pragma once
#ifndef ELASTICROD_H
#define ELASTICROD_H

#include <vector>
#include <Eigen/Dense>
#include "VectorMath.h"
#include "MatrixMath.h"
#include "SceneObjects.h"
#include "SimParameters.h"

// Boundary conodition for the rod
enum BoundaryCondition
{
	BC_FREE,
	BC_FIXED_END,
	BC_RIGIDBODY_END
};

class ElasticRod
{
public:
	ElasticRod();
	ElasticRod(std::vector <Particle, Eigen::aligned_allocator<Particle>> &nodes, SimParameters para, BoundaryCondition bc, std::vector<double> &initialTheta);
	~ElasticRod();

	// Tell the rendering program where a point belongs to a rod element need to be rendered in world coordinate.
	// This program shall interpolate and return the coordinate. The rotation is 
	// the composition of bishop frame of center line and twist at specific location.
	const Eigen::Vector3d renderPos(Eigen::Vector3d &point, int index);


	bool buildConfiguration(Eigen::VectorXd &pos, Eigen::VectorXd &vel);
	bool unbuildConfiguration(const Eigen::VectorXd &pos, const Eigen::VectorXd &vel);

	// For test purpose
	double computeTotalEnergy();
	
	// Velocity Verlet for centerline time integration (unconstrained)
	// And then do step and project
	// This task is separable
	void updateCenterLine();
		
private:
	// Rest position and length
	Eigen::Matrix3Xd restPos;
	Eigen::VectorXd  restLength;
	// Rest material curvature
	Eigen::Matrix2Xd restCurvature;

	// Material Frame (Relatively for 0-th index point)
	Eigen::Vector3d u0, v0, t0;

	// Rod Segments
	std::vector<ElasticRodSegment> rods;

	// Particles
	std::vector<Particle, Eigen::aligned_allocator<Particle>> nodes;

	// Bending Modulus
	// Eigen::Matrix2d B;

	// Twisting Coefficient
	double beta;

	// Use bending stencils to control the binormal curvature and material curvature.
	std::vector <BendingStencil> stencils;

	BoundaryCondition bcStats;

	// Compute theta for each rod segments 
	// Note we know Hessian is tridiagonal, so we use three vectors instead of a sparse matrix
	// since Eigen does not have a tridiagonal linear system solver
	void computeEnergyThetaDifferentialAndHessian(Eigen::VectorXd &dE, Eigen::VectorXd &lowerH, Eigen::VectorXd &centerH, Eigen::VectorXd &upperH);

	// Update rods and stencil information based on new configuration of the position
	void updateAfterTimeIntegration();
	// Update theta
	void updateQuasiStaticFrame();
	// Compute local parallel transport coordinate frame for each rod segments
	void updateBishopFrame();
	// Update all curvature assigned with bending stencils.
	void updateMaterialCurvature();

	SimParameters params;

};


#endif