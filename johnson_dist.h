#ifndef __JOHNSON_DIST_H
#define __JOHNSON_DIST_H

#include "distribution.h"

class JohnsonDistribution : public IDistribution, public IPersistent {
public:
	JohnsonDistribution();
	JohnsonDistribution(double _form, double _shift, double _scale);
	JohnsonDistribution(std::ifstream& file);

	/* ������� ��� ��������� ��������� ����� */
	void setForm(double _form);
	/* ������� ��� ��������� ��������� ������ */
	void setShift(double _shift);
	/* ������� ��� ��������� ��������� �������� */
	void setScale(double _scale);

	/* ������� ��� ��������� ��������� ����� */
	double getForm() const;
	/* ������� ��� ��������� ��������� ������ */
	double getShift() const;
	/* ������� ��� ��������� ��������� �������� */
	double getScale() const;

	/* ��������� ��������� �������� �������������� �� ������ �������� */
	double getRandomVariable() const override;
	/* ���������� ������� ��������� ��� ������������� ��������*/
	double calculateDensity(double x) const override;
	/* ���������� ��������������� �������� ��� ������������� �������� */
	double calculateMathExpectation() const override;
	/* ���������� ��������� ��� ������������� �������� */
	double calculateVariance() const override;
	/* ���������� ����������� �������� ��� ������������� �������� */
	double calculateCoeffKurtosis() const override;
	/* ���������� ����������� ���������� ��� ������������� �������� */
	double calculateCoeffAsymmetry() const override;

	/* ������� ��� ���������� ���������� ������������ �������� � ���� */
	void save(std::ofstream& file) override;
	/* ������� ��� �������� ���������� ������������ �������� �� ����� */
	void load(std::ifstream& file) override;
	void saveDataGraph(const std::vector<double> selection, std::ofstream& file) const override;
	~JohnsonDistribution() {};

private:
	double form;
	double shift;
	double scale;
	/* �������� �������� �� ������������� ����������� */
	bool isStandartDistribution() const;
	/* ��������� ���������� �������������� ��������� �������� �� �������(0; 1) */
	double getUniformRandomVariable() const;
};

#endif // !__JOHNSON_DIST_H
