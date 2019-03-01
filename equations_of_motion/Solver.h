#pragma once

#include <iostream>
#include <complex>
#include <armadillo>
#include <boost/numeric/odeint.hpp>

#include "../auxiliary_functions/double_to_int.h"
#include "../data_types/Tetracube.h"

#include "time_dependent_hamiltonian.h"

using namespace boost::numeric::odeint;
using namespace std::complex_literals;

using state_t = arma::cx_mat;

namespace boost {
    namespace numeric {
        namespace odeint {
            template<>
            struct vector_space_norm_inf<arma::cx_mat> {
                typedef double result_type;
                result_type operator()(const arma::cx_mat &u) const {
                    return arma::norm(u, "inf");
                }
            };
        }}}


template <typename Hamiltonian>
class System{
    using ham_t = Hamiltonian;
public:
    System(const ham_t &h) : h(h) {}
    void operator()(const state_t &u, state_t &du, double t){
            du = -I * h(t) * u;
    }
private:
    const ham_t &h;
};

class Saver{
public:
    Saver(cx_tetracube &u_matrices, double dt) : u_matrices(u_matrices), dt(dt) {}
    int j = 0;
    void operator()(const state_t &u , double t){
//        std::cout << "Saving element (" << t_to_index(t) << ", " << j <<  ")" << std::endl;
//        std::cout << "t_i = " << t << std::endl;
        u_matrices.slice(t_to_index(t), j) = u;
    }
private:
    cx_tetracube &u_matrices;
    const double dt;
    inline int t_to_index(double t){
        return int_round(t / dt);
    }
};

template <typename Hamiltonian>
class Solver{
    using ham_t = Hamiltonian;
    using system_t = System<ham_t>;
    using observer_t = Saver;
    //using stepper_t = runge_kutta4<state_t, double, state_t, double, vector_space_algebra>;
    //using stepper_t = runge_kutta_dopri5<state_t, double, state_t, double, vector_space_algebra>;
    using stepper_t = runge_kutta_cash_karp54<state_t, double, state_t, double, vector_space_algebra>;
    //using stepper_t = adams_bashforth_moulton<5, state_t, double, state_t, double, vector_space_algebra>;
    //using stepper_t = dense_output_runge_kutta<controlled_runge_kutta<runge_kutta_dopri5<state_t, double, state_t, double, vector_space_algebra>>>;
public:
    Solver(int grid_size, double t_max, const ham_t &h, cx_tetracube &u_matrices) : dt(calculate_dt(t_max, grid_size)),
            u_matrices(u_matrices), h(h), equation(h), saver(u_matrices, dt) {}
    void solve_column(int i_0, int i_1, int j, state_t &u_0){
        saver.j = j;
        double dt_with_sign = i_1 < i_0 ? -dt : dt;
        integrate_const(stepper, equation, u_0, index_to_t(i_0), index_to_t(i_1) + dt_with_sign / 10., dt_with_sign, saver);
    }

private:
    const double dt;
    inline double index_to_t(int i) const{
        return i * dt;
    }
    static double calculate_dt(double t_max, int grid_size) {
        std::cout << "solver_step: " << t_max / (grid_size - 1) << std::endl;
        return t_max / (grid_size - 1);
    }
    cx_tetracube &u_matrices;
    const ham_t &h;
    const stepper_t stepper;
    const system_t equation;
    observer_t saver;
};


