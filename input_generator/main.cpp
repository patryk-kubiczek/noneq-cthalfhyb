#include <iostream>
#include <armadillo>
#include "InputGenerator.h"
#include <cmath>

int main(int argc, char **argv)
{
    int n_times = 252;
    double dt = 1. / (250 - 1);


    double beta = 1;
    double t_max = 1;
    int real_grid_size = 250;
    int imag_grid_size = 250;

    assert((n_times - 1) * dt >= t_max);

    int n_flavor = 2;

    std::string name = "input/test";

    int n_bath = 1;
    
    double U = 3;
    double eps_0 = 1;
    double eps = -1;
    double eps_bath = 0;
    double V = 1; //std::sqrt(2 * 1 / (PI * n_bath) * 1);
    
    auto time_dependent_eps_bath = std::vector<double>(n_times, eps_bath);
    auto time_dependent_V = std::vector<cx_double>(n_times, cx_double(V, 0));
    auto time_dependent_eps = std::vector<double>(n_times, eps);
    time_dependent_eps[0] = eps_0;
    auto time_dependent_U = std::vector<double>(n_times, U);

    vec_re_3_t bath_energies;
    vec_3_t hybridizations;
    for(int n = 0; n < n_bath; ++n){
        bath_energies.emplace_back(vec_re_2_t(n_flavor, time_dependent_eps_bath));
        hybridizations.emplace_back(vec_2_t(n_flavor, time_dependent_V));
    }
    vec_re_2_t local_energies(n_flavor, time_dependent_eps);
    vec_re_1_t Hubbard_U(time_dependent_U);

    assert(bath_energies.size() == n_bath);
    assert(hybridizations.size() == n_bath);
    assert(bath_energies[0].size() == n_flavor);
    assert(hybridizations[0].size() == n_flavor);
    assert(local_energies.size() == n_flavor);
    assert(bath_energies[0][0].size() == n_times);
    assert(hybridizations[0][0].size() == n_times);
    assert(local_energies[0].size() == n_times);
    assert(Hubbard_U.size() == n_times);


    UserInput user_input = {
        bath_energies,       // vector in k-space of vectors in flavor space of vectors representing time dependence of bath energies
        hybridizations,      // vector in k-space of vectors in flavor space of vectors representing time dependence of hybridizations
        local_energies,      // vector in flavor space of vectors representing time dependence of local energies
        Hubbard_U,           // vector representing time dependence of Hubbard U
        dt,                  // time resolution of USER input data, data[i] = data(t=i*dt), applies to ALL vectors above
        n_times,             // number of time points in USER input data, applies to ALL vectors above
        beta,                // temperature^{-1}
        t_max,               // maximal time on the real axis
        real_grid_size,      // desired number of real time points in the GENERATED input (i.e. output of this code)
        imag_grid_size,      // desired number of imaginary time points in the GENERATED input (i.e. output of this code)
        name                 // user-defined name of the parameter set, will be used for naming the output
    };

    InputGenerator input_generator = {user_input};

    input_generator.generate_all_for_determinant();
    input_generator.generate_all_for_one_spin();
    input_generator.generate_benchmark_data();
}

