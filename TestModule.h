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

class TestModule
{
public:
	bool testEnergyDifferential(std::function<double()> computeEnergy, std::function<void(Eigen::VectorXd &)> computeGradient);
	bool testEnergyHessian(std::function<void(Eigen::VectorXd &)> computeGradient, std::function<void(Eigen::SparseMatrix<double> &)> computeHessian);
	bool testEnergyHessian(std::function<void(Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> computeGradientandHessian);

};


#endif // !TESTMODULE_H