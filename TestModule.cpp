#include <iostream>

#include "TestModule.h"

bool TestModule::testEnergyDifferential(const Eigen::VectorXd &q, std::function<double(const Eigen::VectorXd &)> computeEnergy, std::function<void(const Eigen::VectorXd &, Eigen::VectorXd &)> computeGradient)
{
	Eigen::VectorXd f;
	computeGradient(q, f);
	std::cout << "Derivative is : " << f.norm() << std::endl;
	double E0;
	E0 = computeEnergy(q);
	std::cout << "Original position is" << q.transpose() << std::endl;
	Eigen::VectorXd direction = Eigen::VectorXd::Random(q.size());
	direction.normalize();

	for (int k = 3; k <= 12; k++)
	{
		double eps = pow(10, -k);
		std::cout << "New proposed : " << (q + eps * direction).transpose() << std::endl;
		double E1 = computeEnergy(q + eps * direction);

		std::cout << "EPS is: " << eps << std::endl;
		std::cout << "Difference of Energy is: " << (E1-E0)/ eps << std::endl;
		std::cout << "Norm of Directional Gradient is: " << -f.dot(direction) << std::endl;
		std::cout << "The difference between above two is: " << (E1 - E0) / eps + f.dot(direction) << std::endl << std::endl;
	}

	return false;
}

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

bool TestModule::testEnergyHessian(const Eigen::VectorXd &q, std::function<void(const Eigen::VectorXd &, Eigen::VectorXd&, Eigen::SparseMatrix<double>&)> computeGradientandHessian)
{
	Eigen::SparseMatrix<double> gradF;
	Eigen::VectorXd f;
	computeGradientandHessian(q, f, gradF);

	for (int k = 3; k <= 12; k++)
	{
		double eps = pow(10, -k);
		Eigen::VectorXd direction = Eigen::VectorXd::Random(q.size());
		direction.normalize();
		Eigen::VectorXd epsF;
		Eigen::SparseMatrix<double>  epsgradF;
		computeGradientandHessian(q + eps * direction, epsF, epsgradF);

		std::cout << "EPS is: " << eps << std::endl;
		std::cout << "Norm of Finite Difference is: " << (epsF - f).norm() / eps << std::endl;
		std::cout << "Norm of Directinal Gradient is: " << (gradF*direction).norm() << std::endl;
		std::cout << "The difference between above two is: " << ((epsF - f) / eps - gradF*direction).norm() << std::endl << std::endl;
	}

	return false;
}
