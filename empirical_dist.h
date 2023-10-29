#ifndef __EMPIRICAL_DIST_H
#define __EMPIRICAL_DIST_H

#include "distribution.h"

class EmpiricalDistribution : public IDistribution, public IPersistent {
public:
	EmpiricalDistribution(int _n, const IDistribution& _d, int _k = 1);
	EmpiricalDistribution(std::ifstream& file);
	EmpiricalDistribution& operator=(const EmpiricalDistribution& d);
	EmpiricalDistribution(const EmpiricalDistribution& d);
	
	/* ������� ��� ��������� ��������� ������� ������� */
	int getN() const;
	/* ������� ��� ��������� ��������� k */
	int getK() const;
	/* ������� ��� ��������� ������� */
	std::vector<double> getSelection() const;
	/* ������� ��� ��������� ������� ���������� */
	std::vector<double> getFrequencies() const;
	/* ������� ��� ��������� ��������� k */

	void setK(int _k);

	/* ��������� ��������� ��������, ������� ������������ ������������� */
	double getRandomVariable() const override;
	/* ���������� ������� ��������� ��� ������������� ������������� */
	double calculateDensity(double x) const override;
	/* ���������� ��������������� �������� ��� ������������� ������������� */
	double calculateMathExpectation() const override;
	/* ���������� ��������� ��� ������������� ��� ������������� ������������� */
	double calculateVariance() const override;
	/* ���������� ����������� ���������� ��� ������������� ������������� */
	double calculateCoeffAsymmetry() const override;
	/* ���������� ����������� �������� ��� ������������� ������������� */
	double calculateCoeffKurtosis() const override;

	/* ������� ��� ���������� ���������� ������������� ������������ � ���� */
	void save(std::ofstream& file) override;
	/* ������� ��� �������� ���������� ������������� ������������ �� ����� */
	void load(std::ifstream& file) override;
	void saveDataGraph(const std::vector<double> selection, std::ofstream& file) const override;
	~EmpiricalDistribution();

private:
	int n;
	int k;
	std::vector<double> selection;
	std::vector<double> boundaries;
	std::vector<double> frequencies;
	/* ���������� k �� ������� ���������� */
	int calculateK() const;
	/* ������������� ������� ��������� ������� */
	std::vector<double> generateSelection(const IDistribution& d);
	/* ���������� ����� ��� ��������� */
	double calculateDelta() const;
	/* ���������� ������� �� ��������� */
	std::vector<double> divideSelectionIntoIntervals() const;
	/* ������������� �������� ����� ��� ������ ������� ��������, �������������� ��������� Xni */
	int leftsideBirnarySearch(double min) const;
	/* �������������� �������� ����� ��� ������ ���������� ��������, �������������� ��������� Xni */
	int rightsideBirnarySearch(double max) const;
	/* ������� ��� ���������� ������� ������������ */
	std::vector<double> calculateFrequency() const;
	/* ����� ���������, �������� ����������� x */
	int getIndexInterval(double x) const;
	/* ���������� ������������ ����������� */
	double calculateCumulProb(int i) const;
};

#endif // !__EMPIRICAL_DIST_H