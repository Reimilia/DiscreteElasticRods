#pragma once
#ifndef ELASTICROD_H
#define ELASTICROD_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "VectorMath.h"
#include "MatrixMath.h"
#include "SceneObjects.h"
#include "SimParameters.h"

// Boundary conodition for the rod


class ElasticRod
{
public:
	ElasticRod();
	ElasticRod(std::vector <Particle, Eigen::aligned_allocator<Particle>> &nodes, SimParameters para);
	~ElasticRod();

	// Assign Boundary Condition
	//bool assignClampedBoundaryCondition(double theta0, double theta1);
	//bool assignRigidBodyBoundaryCondition(int idx0, int idx1);
	bool setSimulationParameters(SimParameters para);

	// Tell the rendering program where a point belongs to a rod element need to be rendered in world coordinate.
	// This program shall interpolate and return the coordinate. The rotation is 
	// the composition of bishop frame of center line and twist at specific location.
	const Eigen::Vector3d renderPos(Eigen::Vector3d &point, int index);


	bool buildConfiguration(Eigen::VectorXd &pos, Eigen::VectorXd &vel);
	bool unbuildConfiguration(const Eigen::VectorXd &pos, const Eigen::VectorXd &vel);

	// For test purpose
	double computeTotalEnergy();

	// Assemble Centerline Forces from potential energy
	void computeCenterlineForces(Eigen::VectorXd &);
		
	// Update rods and stencil information based on new configuration of the position
	void updateAfterPosChange();

	// Compute theta for each rod segments 
	// Note we know Hessian is tridiagonal, so we use three vectors instead of a sparse matrix
	// since Eigen does not have a tridiagonal linear system solver, I use the one in this webpage:
	// https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	void computeEnergyThetaDifferentialAndHessian(const Eigen::VectorXd &theta,Eigen::VectorXd &dE, Eigen::VectorXd &lowerH, Eigen::VectorXd &centerH, Eigen::VectorXd &upperH);

	
	// Rigid body id, if free or clamped only, -1 will be put here;
	int leftRigidBody, rightRigidBody;

	// Template point for two ends
	Eigen::Vector3d leftTemplateCoord, rightTemplateCoord;

	// Rest position and length
	Eigen::MatrixXd restPos;
	Eigen::VectorXd  restLength;
	// Rest material curvature
	Eigen::MatrixXd restCurvature;

	const Particle & getNodeSegment(int index) { return nodes[index]; }
	const ElasticRodSegment & getRodSegment(int index) { return rods[index]; }
	const BendingStencil & getStencilSegment(int index) { return stencils[index]; }

	inline int getNodeNumbers() { return (int)nodes.size(); }
	inline int getRodNumbers() { return (int)rods.size(); }
	inline int getStencilNumbers() { return (int)stencils.size(); }


	
	// Update theta
	bool updateQuasiStaticFrame();
	// Compute local parallel transport coordinate frame for each rod segments
	void updateBishopFrame();
	// Update all curvature assigned with bending stencils.
	void updateMaterialCurvature();

private:
	// Material Frame (Relatively for 0-th index point)
	Eigen::Vector3d u0, v0, t0;

	
	// Particles
	std::vector<Particle, Eigen::aligned_allocator<Particle>> nodes;

	// Rod Segments
	std::vector<ElasticRodSegment> rods;

	// Use bending stencils to control the binormal curvature and material curvature.
	std::vector <BendingStencil> stencils;

	// Twisting Coefficient
	double beta;
	


	SimParameters params;


};


#endif