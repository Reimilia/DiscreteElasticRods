#include "TestModule.h"

bool TestModule::testEnergyDifferential(std::function<double()> computeEnergy, std::function<void(Eigen::VectorXd&)> computeGradient)
{
	return false;
}

bool TestModule::testEnergyHessian(std::function<void(Eigen::VectorXd&)> computeGradient, std::function<void(Eigen::SparseMatrix<double>&)> computeHessian)
{
	return false;
}

bool TestModule::testEnergyHessian(std::function<void(Eigen::VectorXd&, Eigen::SparseMatrix<double>&)> computeGradientandHessian)
{
	return false;
}
