#pragma once

#include <cassert>
#include <armadillo>
#include <vector>
#include <complex>

#include "../model/params.h"
#include "../equations_of_motion/EquationsOfMotion.h"
#include "../auxiliary_functions/save_functions.h"



using vec_1_t = std::vector<cx_double>;
using vec_2_t = std::vector<vec_1_t>;
using vec_3_t = std::vector<vec_2_t>;

using vec_re_1_t = std::vector<double>;
using vec_re_2_t = std::vector<vec_re_1_t>;
using vec_re_3_t = std::vector<vec_re_2_t>;

struct UserInput {
    vec_re_3_t bath_energies;       // vector in k-space of vectors in flavor space of vectors representing time dependence of bath energies
    vec_3_t hybridizations;         // vector in k-space of vectors in flavor space of vectors representing time dependence of hybridizations
    vec_re_2_t local_energies;      // vector in flavor space of vectors representing time dependence of local energies
    vec_re_1_t Hubbard_U;           // vector representing time dependence of Hubbard U
    double dt;                      // time resolution of USER input data, data[i] = data(t=i*dt), applies to ALL vectors above
    int n_times;                    // number of time points in USER input data, applies to ALL vectors above
    double beta;                    // temperature^{-1}
    double t_max;                   // maximal time on the real axis
    int real_grid_size;             // desired number of real time points in the GENERATED input (i.e. output of this code)
    int imag_grid_size;             // desired number of imaginary time points in the GENERATED input (i.e. output of this code)
    std::string params_name;        // user-defined name of the parameter set, will be used for naming the output
};

class InputGenerator {
public:
    InputGenerator(const UserInput &data);
    InputGenerator(vec_re_3_t &bath_energies, 
                   vec_3_t &hybridizations,
                   vec_re_2_t &local_energies,
                   vec_re_1_t &Hubbard_U,
                   double dt,
                   int n_times,
                   double beta,
                   double t_max,
                   int real_grid_size,
                   int imag_grid_size,
                   std::string params_name) : 
                   InputGenerator(UserInput{
                           bath_energies,
                           hybridizations,
                           local_energies,
                           Hubbard_U,
                           dt,
                           n_times,
                           beta,
                           t_max,
                           real_grid_size,
                           imag_grid_size,
                           params_name
                   }) {}

    void generate_p_t_t();
    void generate_p_tau_tau();
    void generate_p() {
        generate_p_t_t();
        generate_p_tau_tau();
    }
    void generate_delta();
    void generate_c_and_a_operator_matrices();
    void generate_all_for_determinant() {
        generate_p();
        generate_delta();
        generate_c_and_a_operator_matrices();
    }

    void generate_u_t_t(int n);
    void generate_u_tau_tau(int n);
    void generate_phi_1_t_t();
    void generate_phi_1_tau_tau();
    void generate_all_for_one_spin(){
        generate_delta();
        generate_u_t_t(0);
        generate_u_t_t(1);
        generate_u_tau_tau(0);
        generate_u_tau_tau(1);
        generate_phi_1_t_t();
        generate_phi_1_tau_tau();
    }
    void generate_noninteracting_benchmark_data();
    void generate_interacting_benchmark_data();
    void generate_benchmark_data(){
        generate_noninteracting_benchmark_data();
        generate_interacting_benchmark_data();
    }

private:
    UserInput data;
};





