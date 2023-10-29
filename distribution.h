#ifndef __DISTRIBUTION_H
#define __DISTRIBUTION_H
#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>

class IDistribution {
public:
	double virtual getRandomVariable() const = 0;
	double virtual calculateDensity(double x) const = 0;
	double virtual calculateMathExpectation() const = 0;
	double virtual calculateVariance() const = 0;
	double virtual calculateCoeffKurtosis() const = 0;
	double virtual calculateCoeffAsymmetry() const = 0;
};

class IPersistent {
	void virtual save(std::ofstream& file) = 0;
	void virtual load(std::ifstream& file) = 0;
	void virtual saveDataGraph(const std::vector<double> selection, std::ofstream& file) const = 0;
};

#endif // !__DISTRIBUTION_H
