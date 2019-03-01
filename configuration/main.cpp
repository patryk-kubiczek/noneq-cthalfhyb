#include "Configuration.h"

int main(){
    double beta = 1;
    double t_max = 1;
    int n_flavor = 2;
    
    RandomNumberGenerator rng;
    rng.initialize(907);
    
    Configuration conf(t_max, beta, n_flavor);
    conf.randomly_generate(3, rng);
    
    std::cout << conf << std::endl;

    int n_ini;
    auto contour_times = conf.get_contour_times(&n_ini);

    std::cout << n_ini << std::endl;
    for(const auto& ct : contour_times){
        std::cout << ct.s(t_max) << std::endl;
    }
}

