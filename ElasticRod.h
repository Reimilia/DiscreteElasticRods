#pragma once
#ifndef ELASTICROD_H
#define ELASTICROD_H

#include <vector>
#include "VectorMath.h"
#include "SceneObjects.h"

class ElasticRod
{
public:
	ElasticRod(int numSegments, Eigen::Vector3d st, Eigen::Vector3d ed);
	~ElasticRod();

	// Build Configuration Dof for time integrator
	void buildconfiguration(Eigen::VectorXd &pos, Eigen::VectorXd &theta, Eigen::VectorXd &vel);
	void updateconfiguration(Eigen::VectorXd &pos, Eigen::VectorXd &theta, Eigen::VectorXd &vel);

private:
	// Material Frame (Relatively for 0-th index point)
	Eigen::Vector3d u0, v0, t0;

	// Rod Segments
	std::vector<ElasticRod> rods;

	// Particles
	std::vector<Particle, Eigen::aligned_allocator<Particle>> nodes;

	// Bending Modulus
	Eigen::Matrix2d B;

	// Twisting Coefficient
	double beta;


	// Boundary conodition for the rod
	enum BoundaryCondition
	{
		BC_FREE,
		BC_FIXED_END,
	};

	// Compute theta for each rod segments 
	void updateQuasiStaticFrame();
	// Compute local parallel transport coordinate frame for each rod segments
	void updateBishopFrame();

};


#endif