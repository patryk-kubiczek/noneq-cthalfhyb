#include "RandomNumberGenerator.h"

void RandomNumberGenerator::initialize(unsigned int seed, unsigned int leap, unsigned int id) {
    this->leap = leap;
    generator.seed(seed);
    if(id > 0) generator.discard(id);
}

double RandomNumberGenerator::random_number(double min, double max) {
    random_number_dist.param(std::uniform_real_distribution<double>::param_type{min, max});
    if(leap > 1) generator.discard(leap - 1);
    return random_number_dist(generator);
}

int RandomNumberGenerator::random_int(int min, int max) {
    random_int_dist.param(std::uniform_int_distribution<int>::param_type{min, max});
    if(leap > 1) generator.discard(leap - 1);
    return random_int_dist(generator);
}

contour_time RandomNumberGenerator::random_contour_time(const contour_time &ct_min, const contour_time &ct_max, double t_max, double beta) {
    double s;
    if(ct_max > ct_min){
        s = random_number(ct_min.s(t_max), ct_max.s(t_max));
    }
    else{
        double contour_length = 2 * t_max + beta;
        s = random_number(ct_min.s(t_max), ct_max.s(t_max) + contour_length);
        if(s > contour_length)
            s -= contour_length;
    }
    assert(s >= 0 && s <= 2 * t_max + beta);
    return contour_time::make_contour_time(s, t_max);
}
