#ifndef __EMPIRICAL_DIST_H
#define __EMPIRICAL_DIST_H

#include "distribution.h"

class EmpiricalDistribution : public IDistribution, public IPersistent {
public:
	EmpiricalDistribution(int _n, const IDistribution& _d, int _k = 1);
	EmpiricalDistribution(std::ifstream& file);
	EmpiricalDistribution& operator=(const EmpiricalDistribution& d);
	EmpiricalDistribution(const EmpiricalDistribution& d);
	
	/* Функция для получения параметра размера выборки */
	int getN() const;
	/* Функция для получения параметра k */
	int getK() const;
	/* Функция для получения выборки */
	std::vector<double> getSelection() const;
	/* Функция для получения массива частностей */
	std::vector<double> getFrequencies() const;
	/* Функция для получения параметра k */

	void setK(int _k);

	/* Генерация случайной величины, имеющей эмпирическое распределение */
	double getRandomVariable() const override;
	/* Вычисление функции плотности для эмпирического распределения */
	double calculateDensity(double x) const override;
	/* Вычисление математического ожидания для эмпирического распределения */
	double calculateMathExpectation() const override;
	/* Вычисление дисперсии для распределения для эмпирического распределения */
	double calculateVariance() const override;
	/* Вычисление коэффицинта асимметрии для эмпирического распределения */
	double calculateCoeffAsymmetry() const override;
	/* Вычисление коэффицинта эксцесса для эмпирического распределения */
	double calculateCoeffKurtosis() const override;

	/* Функция для сохранения параметров эмпирического рапсредления в файл */
	void save(std::ofstream& file) override;
	/* Функция для загрузки параметров эмпирического рапределения из файла */
	void load(std::ifstream& file) override;
	void saveDataGraph(const std::vector<double> selection, std::ofstream& file) const override;
	~EmpiricalDistribution();

private:
	int n;
	int k;
	std::vector<double> selection;
	std::vector<double> boundaries;
	std::vector<double> frequencies;
	/* Вычисление k по формуле Стерджесса */
	int calculateK() const;
	/* Сгенерировать выборку случайных величин */
	std::vector<double> generateSelection(const IDistribution& d);
	/* Вычисление длины для интервала */
	double calculateDelta() const;
	/* Разделение выборки на интервалы */
	std::vector<double> divideSelectionIntoIntervals() const;
	/* Левосторонний бинарный поиск для поиска первого элемента, принадлежащего интервалу Xni */
	int leftsideBirnarySearch(double min) const;
	/* Правосторонний бинарный поиск для поиска последнего элемента, принадлежащего интервалу Xni */
	int rightsideBirnarySearch(double max) const;
	/* Функция для построения массива частотностей */
	std::vector<double> calculateFrequency() const;
	/* Поиск интервала, которому принадлежит x */
	int getIndexInterval(double x) const;
	/* Вычисление кумулятивной вероятности */
	double calculateCumulProb(int i) const;
};

#endif // !__EMPIRICAL_DIST_H