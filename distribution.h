#ifndef __DISTRIBUTION_H
#define __DISTRIBUTION_H
#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>

class IDistribution {
public:
	/* Генерация случайной величины */
	double virtual getRandomVariable() const = 0;
	/* Вычисление функции плотности */
	double virtual calculateDensity(double x) const = 0;
	/* Вычисление математического ожидания */
	double virtual calculateMathExpectation() const = 0;
	/* Вычисление дисперсии */
	double virtual calculateVariance() const = 0;
	/* Вычисление коэффицинта эксцесса */
	double virtual calculateCoeffKurtosis() const = 0;
	/* Вычисление коэффицинта асимметрии */
	double virtual calculateCoeffAsymmetry() const = 0;
};

class IPersistent {
	/* Сохранение распределения в файл */
	void virtual save(std::ofstream& file) = 0;
	/* Загрузка распределения из файла */
	void virtual load(std::ifstream& file) = 0;
	/* Сохранение данных для построения функции плотности */
	void virtual saveDataGraph(const std::vector<double> selection, std::ofstream& file) const = 0;
};

#endif // !__DISTRIBUTION_H
