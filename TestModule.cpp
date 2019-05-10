#include <iostream>

#include "TestModule.h"

bool TestModule::testEnergyDifferential(const Eigen::VectorXd &q, std::function<double(const Eigen::VectorXd &)> computeEnergy, std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &)> computeGradient)
{
	Eigen::VectorXd f;
	computeGradient(q, f);

	double E0;
	E0 = computeEnergy(q);

	for (int k = 1; k <= 12; k++)
	{
		double eps = pow(10, -k);
		Eigen::VectorXd direction = Eigen::VectorXd::Random(q.size());
		direction.normalize();
		double E1 = computeEnergy(q + eps * direction);

		std::cout << "EPS is: " << eps << std::endl;
		std::cout << "Difference of Energy is: " << (E1-E0)/ eps << std::endl;
		std::cout << "Norm of Directional Gradient is: " << f.dot(direction) << std::endl;
		std::cout << "The difference between above two is: " << (E1 - E0) / eps + f.dot(direction) << std::endl << std::endl;
	}
}

bool TestModule::testEnergyHessian(std::function<void(Eigen::VectorXd&)> computeGradient, std::function<void(Eigen::SparseMatrix<double>&)> computeHessian)
{
	return false;
}

bool TestModule::testEnergyHessian(std::function<void(Eigen::VectorXd&, Eigen::SparseMatrix<double>&)> computeGradientandHessian)
{
	return false;
}
