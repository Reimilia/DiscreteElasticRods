#include "ElasticRod.h"


ElasticRod::ElasticRod()
{
}

ElasticRod::ElasticRod(std::vector<ElasticRod> rods, std::vector<Particle, Eigen::aligned_allocator<Particle>> nodes)
{
}

ElasticRod::~ElasticRod()
{
}

const Eigen::Vector3d ElasticRod::renderPos(Eigen::Vector3d & point, int index)
{
}


void ElasticRod::updateQuasiStaticFrame()
{
	/*
	Using Newton's method to solve non-linear equation to determine twist angle theta in each rod element
	Note that it shall incoporates boundary conditions
	*/

	

}

void ElasticRod::updateBishopFrame()
{
	/*	
	Compute Discrete Parallel Transportation to update Bishop Frame for centerline
	Move coordinate frame via the rotational matrix
	*/

	int nRods = (int)rods.size();
	// The first one need to set up manually
	t0 = nodes[1].pos - nodes[0].pos;
	t0 = t0 / t0.norm();
	// u0 $\perp$ t0
	u0 = Eigen::Vector3d(t0[2]-t0[1],t0[0]-t0[2],t0[1]-t0[0]);
	u0 = u0 / u0.norm();
	v0 = u0.cross(t0);

	rods[0].t0 = t0;
	rods[0].u0 = u0;
	rods[0].v0 = v0;


	// Now compute Bishop frame
	for (int i = 1; i < nRods; i++)
	{
		rods[i].t0 = nodes[i + 1].pos - nodes[i].pos;
		rods[i].t0 = (rods[i].t0) / (rods[i].t0).norm();
		Eigen::Vector3d n = (rods[i-1].t0).cross(rods[i].t0);

		rods[i].u0 = VectorMath::rotationMatrix(n*asin(n.norm()) / n.norm()) * rods[i - 1].u0;
		rods[i].u0 = (rods[i].u0) / (rods[i].u0).norm();

		rods[i].v0 = (rods[i].u0).cross(v0);
	}
}
