#pragma once
#ifndef ELASTICROD_H
#define ELASTICROD_H

#include <vector>
#include "VectorMath.h"
#include "MatrixMath.h"
#include "SceneObjects.h"

class ElasticRod
{
public:
	ElasticRod();
	ElasticRod(std::vector<ElasticRod> rods, std::vector < Particle, Eigen::aligned_allocator<Particle>> nodes);
	~ElasticRod();

	// Tell the rendering program where a point belongs to a rod element need to be rendered in world coordinate.
	// This program shall interpolate and return the coordinate. The rotation is 
	// the composition of bishop frame of center line and twist at specific location.
	static const Eigen::Vector3d renderPos(Eigen::Vector3d &point, int index);

private:
	// Rest position and velocity
	Eigen::MatrixX3d restPos;
	Eigen::MatrixX3d restVel;
	// Rest material curvature
	Eigen::MatrixX2d restCurvature;

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


	// Boundary conodition for the rod
	enum BoundaryCondition
	{
		BC_FREE,
		BC_FIXED_END,
		BC_RIGIDBODY_END
	};

	BoundaryCondition bcStats;

	// Compute theta for each rod segments 
	// Note we know Hessian is tridiagonal, so we use three vectors instead of a sparse matrix
	// since Eigen does not have a tridiagonal linear system solver
	void computeEnergyDifferentialAndHessian(Eigen::VectorXd &dE, Eigen::VectorXd &lowerH, Eigen::VectorXd &centerH, Eigen::VectorXd &upperH);
	void updateQuasiStaticFrame();
	// Compute local parallel transport coordinate frame for each rod segments
	void updateBishopFrame();

};


#endif