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
	static bool testEnergyDifferential(const Eigen::VectorXd &, std::function<double(const Eigen::VectorXd &)> computeEnergy, std::function<void(const Eigen::VectorXd &, Eigen::VectorXd&)> computeGradient);
	static bool testEnergyDifferential(std::function<double()> computeEnergy, std::function<void(Eigen::VectorXd&)> computeGradient);
	static bool testEnergyHessian(std::function<void(Eigen::VectorXd &)> computeGradient, std::function<void(Eigen::SparseMatrix<double> &)> computeHessian);
	static bool testEnergyHessian(std::function<void(Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> computeGradientandHessian);
	static bool testEnergyHessian(const Eigen::VectorXd &, std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &, Eigen::SparseMatrix<double> &)> computeGradientandHessian);

};


#endif // !TESTMODULE_H