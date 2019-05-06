/*
	Test Function that used to check correctness

	Currently two major components are:

	Energy-Force Differential Checker
	Force-Hessian Differential Checker

*/
#ifndef TESTMODULE_H
#define TESTMODULE_H

#include <functional>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace TestModule
{
	bool testEnergyDifferential(std::function<double()> computeEnergy, std::function<void(Eigen::VectorXd &)> computeGradient)
	{
		return false;
	}
	bool testEnergyHessian(std::function<void(Eigen::VectorXd &)> computeGradient, std::function<void(Eigen::SparseMatrix<double> &)> computeHessian)
	{
		return false;
	}
	bool testEnergyHessian(std::function<void(Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> computeGradientandHessian)
	{
		return false;
	}

}


#endif // !TESTMODULE_H