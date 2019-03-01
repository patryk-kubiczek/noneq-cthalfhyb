#pragma once

#include <random>

#include "../model/contour_time.h"

class RandomNumberGenerator {
public:
    void initialize(unsigned int seed, unsigned int leap = 0, unsigned int id = 0);
    double random_number(double min, double max);
    int random_int(int min, int max);
    contour_time random_contour_time(const contour_time &ct_min, const contour_time &ct_max, double t_max, double beta);

private:
    unsigned int leap = 1;
    std::mt19937 generator;
    std::uniform_real_distribution<double> random_number_dist;
    std::uniform_int_distribution<int> random_int_dist;
};