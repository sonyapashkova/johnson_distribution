#include "empirical_dist.h"

int EmpiricalDistribution::calculateK() const {
	return (int)ceil(log2(n) + 1);
}

std::vector<double> EmpiricalDistribution::generateSelection(const IDistribution& d) {
	std::vector<double> selection;
	for (int i = 0; i < n; i++) {
		selection.push_back(d.getRandomVariable());
	}
	sort(selection.begin(), selection.end());
	return selection;
}

double EmpiricalDistribution::calculateDelta() const {
	return (1.0 / k) * (selection[n - 1] - selection[0]);
}

std::vector<double> EmpiricalDistribution::divideSelectionIntoIntervals() const {
	std::vector<double> boundaries;
	double next = selection[0];
	double max = selection[n - 1];
	double delta = calculateDelta();
	boundaries.push_back(next);
	while (next < max) {
		boundaries.push_back(next + delta);
		next += delta;
	}
	return boundaries;
}

int EmpiricalDistribution::leftsideBirnarySearch(double min) const {
	int left = 0;
	int right = n - 1;
	while (left + 1 < right) {
		int midle = left + (right - left) / 2;
		if (selection[midle] >= min) {
			right = midle;
		}
		else {
			left = midle;
		}
	}
	if (selection[right] >= min) {
		return right;
	}
	if (selection[left] >= min) {
		return left;
	}
	return -1;
}

int EmpiricalDistribution::rightsideBirnarySearch(double max) const {
	int left = 0;
	int right = n - 1;
	while (left + 1 < right) {
		int midle = left + floor((right - left) / 2);
		if (selection[midle] < max) {
			left = midle;
		}
		else {
			right = midle;
		}
	}
	if (selection[right] < max) {
		return right;
	}
	if (selection[left] < max) {
		return left;
	}
	return -1;
}

std::vector<double> EmpiricalDistribution::calculateFrequency() const {
	std::vector<double> frequencies;
	int left, right;
	int ni;
	double delta = calculateDelta();
	for (int i = 0; i < boundaries.size() - 1; i++) {
		if (i == 0) {
			right = rightsideBirnarySearch(boundaries[i + 1]);
			ni = right + 1;
		}
		else if (i == boundaries.size() - 2) {
			left = leftsideBirnarySearch(boundaries[i]);
			ni = n - left;
		}
		else {
			left = leftsideBirnarySearch(boundaries[i]);
			right = rightsideBirnarySearch(boundaries[i + 1]);
			ni = right - left + 1;
		}
		frequencies.push_back(ni / (n * delta));
	}
	return frequencies;
}

EmpiricalDistribution::EmpiricalDistribution(int _n, const IDistribution& _d, int _k) :
	n(_n > 1 ? _n : throw 1), k(_k > 1 ? _k : calculateK()), selection(generateSelection(_d)), boundaries(divideSelectionIntoIntervals()), frequencies(calculateFrequency()) {}


EmpiricalDistribution::EmpiricalDistribution(const EmpiricalDistribution& d) {
	n = d.n;
	k = d.k;
	for (int i = 0; i < d.selection.size(); i++) {
		selection.push_back(d.selection[i]);
	}
	for (int i = 0; i < d.boundaries.size(); i++) {
		boundaries.push_back(d.boundaries[i]);
	}
	for (int i = 0; i < d.frequencies.size(); i++) {
		frequencies.push_back(d.frequencies[i]);
	}
}

EmpiricalDistribution& EmpiricalDistribution:: operator=(const EmpiricalDistribution& d) {
	if (this == &d) return *this;
	selection.clear();
	frequencies.clear();
	boundaries.clear();
	n = d.n;
	k = d.k;
	for (int i = 0; i < d.selection.size(); i++) {
		selection.push_back(d.selection[i]);
	}
	for (int i = 0; i < d.boundaries.size(); i++) {
		boundaries.push_back(d.boundaries[i]);
	}
	for (int i = 0; i < d.frequencies.size(); i++) {
		frequencies.push_back(d.frequencies[i]);
	}
	return *this;
}

int EmpiricalDistribution::getN() const {
	return n;
}

int EmpiricalDistribution::getK() const {
	return k;
}

std::vector<double> EmpiricalDistribution::getSelection() const {
	return selection;
}

std::vector<double> EmpiricalDistribution::getFrequencies() const {
	return frequencies;
}

void EmpiricalDistribution::setK(int _k) {
	if (_k <= 1) {
		_k = calculateK();
	}
	k = _k;
	boundaries = divideSelectionIntoIntervals();
	frequencies = calculateFrequency();
}

EmpiricalDistribution::EmpiricalDistribution(std::ifstream& file) {
	if (!file.is_open()) {
		throw 0;
	}
	int _n, _k;
	file >> _n;
	if (_n <= 1) {
		throw 1;
	}
	n = _n;
	for (int i = 0; i < n; i++) {
		double temp;
		file >> temp;
		selection.push_back(temp);
	}
	file >> _k;
	if (_k <= 1) {
		_k = calculateK();
	}
	k = _k;
	boundaries = divideSelectionIntoIntervals();
	frequencies = calculateFrequency();
}

int EmpiricalDistribution::getIndexInterval(double x) const {
	for (int i = 0; i < boundaries.size() - 1; i++) {
		if (x >= boundaries[i] && x < boundaries[i + 1]) {
			return i;
		}
	}
	if (x >= boundaries[boundaries.size() - 2] && x <= boundaries[boundaries.size() - 1]) {
		return boundaries.size() - 2;
	}
}

double EmpiricalDistribution::calculateCumulProb(int k) const {
	double q = 0;
	for (int i = 0; i <= k; ++i) {
		q += frequencies[i];
	}
	return q;
}

double EmpiricalDistribution::getRandomVariable() const {
	double r;
	double topBound = calculateCumulProb(frequencies.size() - 1);
	do r = (double)rand() / RAND_MAX * topBound;
	while (r == 0 || r == topBound);
	for (int i = 0; i < boundaries.size() - 1; ++i) {
		if (r > calculateCumulProb(i) and r < calculateCumulProb(i + 1)) {
			do r = (double)rand() / RAND_MAX * (boundaries[i + 1] - boundaries[i]) + boundaries[i + 1];
			while (r == boundaries[i] || r == boundaries[i + 1]);
			break;
		}
	}
	return r;
}

double EmpiricalDistribution::calculateDensity(double x) const {
	return frequencies[getIndexInterval(x)];
}

double EmpiricalDistribution::calculateMathExpectation() const {
	double sum = 0;
	for (auto& i : selection) {
		sum += i;
	}
	return sum / n;
}

double EmpiricalDistribution::calculateVariance() const {
	double M = calculateMathExpectation();
	double sum = 0;
	for (auto& i : selection) {
		sum += pow(i - M, 2);
	}
	return sum / n;
}

double EmpiricalDistribution::calculateCoeffAsymmetry() const {
	double M = calculateMathExpectation();
	double D = calculateVariance();
	double sum = 0;
	for (auto& i : selection) {
		sum += pow(i - M, 3);
	}
	return sum / (n * pow(D, 1.5));
}

double EmpiricalDistribution::calculateCoeffKurtosis() const {
	double M = calculateMathExpectation();
	double D = calculateVariance();
	double sum = 0;
	for (auto& i : selection) {
		sum += pow(i - M, 4);
	}
	return sum / (n * pow(D, 2)) - 3;
}

void EmpiricalDistribution::save(std::ofstream& file) {
	file << n << "\n";
	for (int i = 0; i < n; i++) {
		file << selection[i] << "\n";
	}
	file << k << "\n";
}

void EmpiricalDistribution::load(std::ifstream& file) {
	selection.clear();
	boundaries.clear();
	frequencies.clear();
	if (!file.is_open()) {
		throw 0;
	}
	int _n, _k;
	file >> _n;
	if (_n <= 1) {
		throw 1;
	}
	n = _n;
	for (int i = 0; i < n; i++) {
		double temp;
		file >> temp;
		selection.push_back(temp);
	}
	file >> _k;
	if (_k <= 1) {
		_k = calculateK();
	}
	k = _k;
	boundaries = divideSelectionIntoIntervals();
	frequencies = calculateFrequency();
}

EmpiricalDistribution::~EmpiricalDistribution() {
	selection.clear();
	boundaries.clear();
	frequencies.clear();
}

void EmpiricalDistribution::saveDataGraph(const std::vector<double> selection, std::ofstream& file) const {
	if (!file.is_open()) {
		throw 0;
	}
	for (int i = 0; i < selection.size(); i++) {
		file << selection[i] << " " << calculateDensity(selection[i]) << "\n";
	}
}