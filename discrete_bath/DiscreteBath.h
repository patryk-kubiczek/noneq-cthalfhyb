#pragma once

#include <vector>
#include <armadillo>
#include <cassert>
#include "../auxiliary_functions/constants.h"
#include "../model/contour_time.h"

class DiscreteBath{
public:
    DiscreteBath( double eps_0, double eps, double U_0, double U,
                  const std::vector<cx_double> &V_0,  const std::vector<cx_double> &V,
                  const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath,
                  double beta, double t_max, int n_blocks);

    void generate_u(const std::vector<contour_time> &time_points, int n_ini);
    void generate_QDT();
    cx_double log_z();
    cx_double log_z_without_QDT();
    arma::cx_mat density_matrix_at_t_max() const;
    arma::cx_mat density_matrix_at_t_max_without_QDT() const;
    
    cx_long_double calculate_z(const std::vector<contour_time> &time_points, int n_ini) {
        generate_u(time_points, n_ini);
        if(n_blocks > 1){
            generate_QDT();
            return std::exp(cx_long_double{log_z()});
        }
        else
            return std::exp(cx_long_double{log_z_without_QDT()});
    }

    arma::cx_mat calculate_density_matrix_at_t_max(int n_ini, const std::vector<contour_time> &time_points, int *n_down_at_max) {
        generate_u(time_points, n_ini);
        if(n_down_at_max) *n_down_at_max = n_at_t_max;
        if(n_blocks > 1) {
            generate_QDT();
            return density_matrix_at_t_max();
        }
        else
            return density_matrix_at_t_max_without_QDT();
    }


private:
    double eps_0 = 0;
    double eps = 0;
    double U_0 = 0;
    double U = 0;

    double beta = 1;
    double t_max = 1;

    int N;
    int n_blocks = 1;

    arma::cx_mat O_R_0;
    arma::cx_mat O_R_1;

    arma::mat O_I_0;
    arma::mat O_I_1;

    arma::cx_mat O_R_01;
    arma::mat O_I_01;
    arma::cx_mat O_IR_0;
    arma::cx_mat O_IR_1;

    arma::vec e_R_0;
    arma::vec e_R_1;
    arma::vec e_I_0;
    arma::vec e_I_1;

    arma::cx_cube u_t;
    arma::cube u_tau;

    arma::cube Q;
    arma::cube R;
    arma::mat D;
    arma::cube T;

    arma::vec D_b;
    arma::vec D_s;

    int n_at_t_max;

    void generate_e_and_O(const std::vector<cx_double> &V_0,  const std::vector<cx_double> &V,
                          const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath);
};
