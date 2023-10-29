#include "catch.hpp"
#include "johnson_dist.h"
#include "empirical_dist.h"
#include "mixture_dist.cpp"


TEST_CASE("[Johnson Distribution] Standart Distribution") {
    JohnsonDistribution d = JohnsonDistribution(1.5, 0, 1);
    CHECK(round(d.calculateVariance() * 1000) / 1000 == 0.716);
    CHECK(round(d.calculateCoeffKurtosis() * 1000) / 1000 == 3.891);
    CHECK(round(d.calculateDensity(0) * 1000) / 1000 == 0.598);
}


TEST_CASE("[Johnson Distribution] Shift-Scale Transformations") {
    JohnsonDistribution d = JohnsonDistribution(2.5, 1, 2);
    CHECK(round(d.calculateVariance() * 1000) / 1000 == 0.754);
    CHECK(round(d.calculateCoeffKurtosis() * 1000) / 1000 == 0.825);
    CHECK(round(d.calculateDensity(0) * 1000) / 1000 == 0.216);
    CHECK(d.calculateMathExpectation() == 1);
}

TEST_CASE("[Mixture Distribution] Trivial Case") {
    JohnsonDistribution d1 = JohnsonDistribution(2.5, 1, 2);
    JohnsonDistribution d2 = JohnsonDistribution(2.5, 1, 2);
    double p = 0.7;
    MixtureDistribution<JohnsonDistribution, JohnsonDistribution> d(d1, d2, p);
    CHECK(round(d.calculateDensity(2) * 1000) / 1000 == 0.216);
    CHECK(d.calculateMathExpectation() == 1);
    CHECK(round(d.calculateVariance() * 1000) / 1000 == 0.754);
    CHECK(d.calculateCoeffAsymmetry() == 0);
    CHECK(round(d.calculateCoeffKurtosis() * 1000) / 1000 == 0.825);
}

TEST_CASE("[Mixture Distribution] Non-Trivial Case") {
    JohnsonDistribution d1 = JohnsonDistribution(2.5, 1, 2);
    JohnsonDistribution d2 = JohnsonDistribution(3, 2, 3);
    double p = 0.5;
    MixtureDistribution<JohnsonDistribution, JohnsonDistribution> d(d1, d2, p);
    CHECK(round(d.calculateDensity(2) * 1000) / 1000 == 0.308);
    CHECK(d.calculateMathExpectation() == 1.5);
    CHECK(round(d.calculateVariance() * 1000) / 1000 == 1.187);
    CHECK(round(d.calculateCoeffAsymmetry() * 1000) / 1000 == 0.212);
    CHECK(round(d.calculateCoeffKurtosis() * 1000) / 1000 == 0.384);
}