#include "johnson_dist.h"

JohnsonDistribution::JohnsonDistribution() :
	form(1.0), shift(0.0), scale(1.0) {}

JohnsonDistribution::JohnsonDistribution(double _form, double _shift, double _scale) :
	form(_form > 0 ? _form : throw 1), shift(_shift), scale(_scale > 0 ? _scale : throw 1) {}


JohnsonDistribution::JohnsonDistribution(std::ifstream& file) {
	double _form, _shift, _scale;
	file.open("johnson.txt");
	if (!file.is_open()) {
		throw 0;
	}
	file >> _form >> _shift >> _scale;
	file.close();
	JohnsonDistribution(_form, _shift, _scale);
}

void JohnsonDistribution::setForm(double _form) {
	if (_form <= 0) {
		throw 1;
	}
	form = _form;
}

void JohnsonDistribution::setShift(double _shift) {
	shift = _shift;
}

void JohnsonDistribution::setScale(double _scale) {
	if (_scale <= 0) {
		throw 1;
	}
	scale = _scale;
}

double JohnsonDistribution::getForm() const {
	return form;
}

double JohnsonDistribution::getShift() const {
	return shift;
}

double JohnsonDistribution::getScale() const {
	return scale;
}

bool JohnsonDistribution::isStandartDistribution() const {
	return shift == 0 && scale == 1;
}

double JohnsonDistribution::getUniformRandomVariable() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

double JohnsonDistribution::getRandomVariable() const {
	double r1 = getUniformRandomVariable();
	double r2 = getUniformRandomVariable();
	double z = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
	double x1 = sinh(z / form);

	if (isStandartDistribution()) {
		return x1;
	}
	return shift + scale * x1;
}

double JohnsonDistribution::calculateDensity(double x) const {
	double x1 = (x - shift) / scale;
	return (1.0 / scale) * (form / sqrt(2 * M_PI)) * (1.0 / sqrt(pow(x1, 2) + 1)) * exp(-pow(form, 2) / 2 * pow(log(x1 + sqrt(pow(x1, 2) + 1)), 2));
}

double JohnsonDistribution::calculateMathExpectation() const {
	if (isStandartDistribution()) {
		return 0;
	}
	return shift;
}

double JohnsonDistribution::calculateVariance() const {
	double w = exp(2 / pow(form, 2));
	double sigma = (w - 1) / 2.0;

	if (isStandartDistribution()) {
		return sigma;
	}
	return pow(scale, 2) * sigma;
}

double JohnsonDistribution::calculateCoeffKurtosis() const {
	double w = exp(2 / pow(form, 2));
	return (pow(w, 2) + 2 * w - 3) / 2.0;
}

double JohnsonDistribution::calculateCoeffAsymmetry() const {
	return 0.0;
}

void JohnsonDistribution::save(std::ofstream& file) {
	file << form << "\n" << shift << "\n" << scale << "\n";
}

void JohnsonDistribution::load(std::ifstream& file) {
	double _form, _shift, _scale;
	if (!file.is_open()) {
		throw 0;
	}
	file >> _form >> _shift >> _scale;
	if (_form <= 0 || _scale <= 0) {
		throw 1;
	}
	form = _form;
	shift = _shift;
	scale = _scale;
}

void JohnsonDistribution::saveDataGraph(const std::vector<double> selection, std::ofstream& file) const {
	if (!file.is_open()) {
		throw 0;
	}
	for (int i = 0; i < selection.size(); i++) {
		file << selection[i] << " " << calculateDensity(selection[i]) << "\n";
	}
}