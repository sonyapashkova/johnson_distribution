#include "distribution.h"

template<class Distribution1, class Distribution2>
class MixtureDistribution : public IDistribution, public IPersistent {
public:
	MixtureDistribution(Distribution1& _d1, Distribution2& _d2, double _p) :
		d1(_d1), d2(_d2), p(_p) {};
	MixtureDistribution(std::ifstream& file);

	Distribution1& component1() { return d1; }
	Distribution2& component2() { return d2; }

	/* Функция для установки параметра смеси */
	void setP(double _p);
	/* Функция для получения параметра смеси */
	double getP() const;

	/* Генерация случайной величины, имеющей распределение в виде смеси двух распределений */
	double getRandomVariable() const override;
	/* Вычисление функции плотности для распределения смесей */
	double calculateDensity(double x) const override;
	/* Вычисление математического ожидания для распределения смесей */
	double calculateMathExpectation() const override;
	/* Вычисление дисперсии для распределения смесей */
	double calculateVariance() const override;
	/* Вычисление коэффицинта асимметрии для распределения смесей */
	double calculateCoeffAsymmetry() const override;
	/* Вычисление коэффицинта эксцесса для распределения смесей */
	double calculateCoeffKurtosis() const override;

	/* Функция для сохранения параметров рапсредления смесей в файл */
	void save(std::ofstream& file) override;
	/* Функция для загрузки параметров рапределения смесей из файла */
	void load(std::ifstream& file) override;
	/* Функция для сохранения данных в файл для построения графика плостности распределения смесей */
	void saveDataGraph(const std::vector<double> selection, std::ofstream& file) const override;
	~MixtureDistribution() {}

private:
	double p;
	Distribution1 d1;
	Distribution2 d2;
	/* Генерация равномерно распределенной случайной величины на отрезке(0; 1) */
	double getUniformRandomVariable() const;
};

template<class dist1, class dist2>
MixtureDistribution<dist1, dist2>::MixtureDistribution(std::ifstream& file) {
	double _p;
	component1().load(file);
	component2().load(file);
	file >> _p;
	if (_p < 0 || _p > 1) {
		throw 1;
	}
	p = _p;
}

template<class dist1, class dist2>
void MixtureDistribution<dist1, dist2>::setP(double _p) {
	if (_p < 0 || _p > 1) {
		throw 1;
	}
	p = _p;
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::getP() const {
	return p;
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::getUniformRandomVariable() const {
	double r;
	do r = (double)rand() / RAND_MAX; while (r == 0 || r == 1);
	return r;
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::getRandomVariable() const {
	double r = getUniformRandomVariable();
	if (r > p) {
		return d1.getRandomVariable();
	}
	return d2.getRandomVariable();
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::calculateDensity(double x) const {
	return (1 - p) * d1.calculateDensity(x) + p * d2.calculateDensity(x);
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::calculateMathExpectation() const {
	return (1 - p) * d1.calculateMathExpectation() + p * d2.calculateMathExpectation();;
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::calculateVariance() const {
	double M = calculateMathExpectation();
	double M1 = d1.calculateMathExpectation();
	double M2 = d2.calculateMathExpectation();
	double D1 = d1.calculateVariance();
	double D2 = d2.calculateVariance();
	return ((1 - p) * (pow(M1, 2) + D1) + p * (pow(M2, 2) + D2)) - pow(M, 2);
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::calculateCoeffAsymmetry() const {
	double M = calculateMathExpectation();
	double M1 = d1.calculateMathExpectation();
	double M2 = d2.calculateMathExpectation();
	double D = calculateVariance();
	double D1 = d1.calculateVariance();
	double D2 = d2.calculateVariance();
	double gamma11 = d1.calculateCoeffAsymmetry();
	double gamma12 = d2.calculateCoeffAsymmetry();
	return (1.0 / pow(D, 1.5)) * ((1 - p) * (pow(M1 - M, 3) + 3 * (M1 - M) * D1 + pow(D1, 1.5) * gamma11) +
		p * (pow(M2 - M, 3) + 3 * (M2 - M) * D2 + pow(D2, 1.5) * gamma12));
}

template<class dist1, class dist2>
double MixtureDistribution<dist1, dist2>::calculateCoeffKurtosis() const {
	double M = calculateMathExpectation();
	double M1 = d1.calculateMathExpectation();
	double M2 = d2.calculateMathExpectation();
	double D = calculateVariance();
	double D1 = d1.calculateVariance();
	double D2 = d2.calculateVariance();
	double gamma11 = d1.calculateCoeffAsymmetry();
	double gamma12 = d2.calculateCoeffAsymmetry();
	double gamma21 = d1.calculateCoeffKurtosis();
	double gamma22 = d2.calculateCoeffKurtosis();
	return (1.0 / pow(D, 2)) * ((1 - p) * (pow(M1 - M, 4) + 6 * pow(M1 - M, 2) * D1 + 4 * (M1 - M) * pow(D1, 1.5) * gamma11 + pow(D1, 2) * (gamma21 + 3)) +
		p * (pow(M2 - M, 4) + 6 * pow(M2 - M, 2) * D2 + 4 * (M2 - M) * pow(D2, 1.5) * gamma12 + pow(D2, 2) * (gamma22 + 3))) - 3;
}

template<class dist1, class dist2>
void MixtureDistribution<dist1, dist2>::save(std::ofstream& file) {
	component1().save(file);
	component2().save(file);
	file << p << "\n";
}

template<class dist1, class dist2>
void MixtureDistribution<dist1, dist2>::load(std::ifstream& file) {
	double _p;
	component1().load(file);
	component2().load(file);
	file >> _p;
	if (_p < 0 || _p > 1) {
		throw 1;
	}
	p = _p;
}

template<class dist1, class dist2>
void MixtureDistribution<dist1, dist2>::saveDataGraph(const std::vector<double> selection, std::ofstream& file) const {
	if (!file.is_open()) {
		throw 0;
	}
	for (int i = 0; i < selection.size(); i++) {
		file << selection[i] << " " << calculateDensity(selection[i]) << "\n";
	}
}