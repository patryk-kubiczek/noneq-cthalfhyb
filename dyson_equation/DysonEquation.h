#pragma once

#include <cassert>
#include <armadillo>
#include <vector>
#include "../auxiliary_functions/constants.h"
#include "../model/contour_time.h"


class DysonEquation{
public:
    DysonEquation(double t_max, double beta, int N_t, int N_tau,
                  double eps_0, double eps, double U_0, double U,
                  const arma::cx_mat &Delta, const arma::cx_vec &Delta_jump);

    DysonEquation(double t_max, double beta, int N_t, int N_tau,
                  double eps_0, double eps, double U_0, double U,
                  const std::vector<cx_double> &V_0,  const std::vector<cx_double> &V,
                  const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath);

    void generate_du(const std::vector<contour_time>& time_points, int n_ini);
    void generate_G0();
    void generate_F_and_M();

    cx_double Z();
    cx_double Z_from_G();
    
    void solve_for_G();
    arma::cx_mat get_G() { return G; }
    arma::cx_vec get_G_jump() { return G0_jump; }
    
    cx_long_double calculate_trace(const std::vector<contour_time>& time_points, int n_ini){
        generate_du(time_points, n_ini);
        void generate_G0();
        void generate_F_and_M();
        return Z();
    }

    cx_double calculate_n_at_t_max(int& n_down_at_t_max){
        n_down_at_t_max = this->n_down_at_t_max;
        solve_for_G();
        return get_G()(i_t_max, i_t_max) + get_G_jump()(i_t_max);
    }
    
    arma::cx_mat calculate_G(int& n_down_at_t_max){
        n_down_at_t_max = this->n_down_at_t_max;
        solve_for_G();
        return get_G() + arma::diagmat(get_G_jump());
    }
    

private:
    const double dt;
    const double dtau;
    const int N_t;
    const int N_tau;
    int N;
    int i_t_max;
    int i_0;

    double eps_0 = 0;
    double eps = 0;
    double U_0 = 0;
    double U = 0;

    std::vector<cx_double> du;

    arma::cx_mat G0;
    arma::cx_vec G0_jump;
    arma::cx_mat Delta;
    arma::cx_vec Delta_jump;

    arma::cx_vec W;

    arma::cx_mat F;
    arma::cx_mat M;
    arma::cx_mat G;

    cx_double Z_bath = 1.;
    cx_double Z0;
    
    int n_down_at_t_max;

    void initialize();
    void generate_w();
    void generate_simple_Delta(const std::vector<cx_double> &V_0,  const std::vector<cx_double> &V,
                               const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath);

    arma::cx_mat convolution(const arma::cx_mat& A, const arma::cx_mat &B, const arma::cx_vec& DA, const arma::cx_vec& DB);



};

