#ifndef __JOHNSON_DIST_H
#define __JOHNSON_DIST_H

#include "distribution.h"

class JohnsonDistribution : public IDistribution, public IPersistent {
public:
	JohnsonDistribution();
	JohnsonDistribution(double _form, double _shift, double _scale);
	JohnsonDistribution(std::ifstream& file);

	/* Функция для установки параметра формы */
	void setForm(double _form);
	/* Функция для установки параметра сдвига */
	void setShift(double _shift);
	/* Функция для установки параметра масштаба */
	void setScale(double _scale);

	/* Функция для получения параметра формы */
	double getForm() const;
	/* Функция для получения параметра сдвига */
	double getShift() const;
	/* Функция для получения параметра масштаба */
	double getScale() const;

	/* Генерация случайной величины распределенной по закону Джонсона */
	double getRandomVariable() const override;
	/* Вычисление функции плотности для распределения Джонсона*/
	double calculateDensity(double x) const override;
	/* Вычисление математического ожидания для распределения Джонсона */
	double calculateMathExpectation() const override;
	/* Вычисление дисперсии для распределения Джонсона */
	double calculateVariance() const override;
	/* Вычисление коэффицинта эксцесса для распределения Джонсона */
	double calculateCoeffKurtosis() const override;
	/* Вычисление коэффицинта асимметрии для распределения Джонсона */
	double calculateCoeffAsymmetry() const override;

	/* Функция для сохранения параметров рапсредления Джонсона в файл */
	void save(std::ofstream& file) override;
	/* Функция для загрузки параметров рапределения Джонсона из файла */
	void load(std::ifstream& file) override;
	void saveDataGraph(const std::vector<double> selection, std::ofstream& file) const override;
	~JohnsonDistribution() {};

private:
	double form;
	double shift;
	double scale;
	/* Проверка является ли распределение стандартным */
	bool isStandartDistribution() const;
	/* Генерация равномерно распределенной случайной величины на отрезке(0; 1) */
	double getUniformRandomVariable() const;
};

#endif // !__JOHNSON_DIST_H
