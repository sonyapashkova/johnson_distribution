#define CATCH_CONFIG_RUNNER

#include <fstream>
#include "catch.hpp"
#include "johnson_dist.h"
#include "empirical_dist.h"
#include "mixture_dist.cpp"

int main(int argc, char* argv[]) {
	setlocale(LC_ALL, "ru");

	JohnsonDistribution d1(4, 0, 1);
	/*JohnsonDistribution d2(3, -3, 1);
	JohnsonDistribution d3(2, 0, 5);
	JohnsonDistribution d4(4, 6, 6);
	MixtureDistribution<JohnsonDistribution, JohnsonDistribution> md1(d1, d2, 0.5);
	MixtureDistribution<JohnsonDistribution, JohnsonDistribution> md2(d3, d4, 0.5);
	MixtureDistribution<JohnsonDistribution, JohnsonDistribution> md3(d1, d2, 0.5);*/
	EmpiricalDistribution ed1(10000, d1);
	EmpiricalDistribution ed2(10000, ed1);
	std::vector<double> selection1 = ed1.getSelection();
	std::vector<double> selection2 = ed2.getSelection();
	std::ofstream file_johnson("johnson_graph.txt");
	std::ofstream file_empirical("empirical_graph.txt");
	std::ofstream file_mixture("mixture_graph.txt");
	d1.saveDataGraph(selection1, file_johnson);
	ed2.saveDataGraph(selection2, file_empirical);
	file_johnson.close();
	file_empirical.close();

	std::cout << "M = " << d1.calculateMathExpectation() << "\n";
	std::cout << "M* = " << ed1.calculateMathExpectation() << "\n";
	std::cout << "M** = " << ed2.calculateMathExpectation() << "\n";
	std::cout << "D = " << d1.calculateVariance() << "\n";
	std::cout << "D* = " << ed1.calculateVariance() << "\n";
	std::cout << "D** = " << ed2.calculateVariance() << "\n";
	std::cout << "gamma1 = " << d1.calculateCoeffAsymmetry() << "\n";
	std::cout << "gamma1* = " << ed1.calculateCoeffAsymmetry() << "\n";
	std::cout << "gamma1** = " << ed2.calculateCoeffAsymmetry() << "\n";
	std::cout << "gamma2 = " << d1.calculateCoeffKurtosis() << "\n";
	std::cout << "gamma2* = " << ed1.calculateCoeffKurtosis() << "\n";
	std::cout << "gamma2** = " << ed2.calculateCoeffKurtosis() << "\n";

	int result = Catch::Session().run(argc, argv);
	return result;
}