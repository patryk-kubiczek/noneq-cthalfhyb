#pragma once

#include "../auxiliary_functions/constants.h"

// Params

struct InputParams{
    double dt = 0.;
    int real_grid_size = 0;                 // number of real time points in the input data
    double dtau = 0.;
    int imag_grid_size = 0;                 // number of imaginary time points in the input data
};

struct ModelParams{
    double eps_0_up;
    double eps_up;
    double eps_0_down;
    double eps_down;

    double U_0;
    double U;

    std::vector<cx_double> V_0_up;
    std::vector<cx_double> V_up;
    std::vector<double> eps_bath_0_up;
    std::vector<double> eps_bath_up;
    std::vector<cx_double> V_0_down;
    std::vector<cx_double> V_down;
    std::vector<double> eps_bath_0_down;
    std::vector<double> eps_bath_down;
};

#include <string>

const std::string extension = "arma";

const std::string p_t_t_name = "p_t_t";
const std::string p_tau_tau_name = "p_tau_tau";
const std::string delta_greater_t_t_name = "delta_greater_t_t";
const std::string delta_lesser_t_t_name = "delta_lesser_t_t";
//const std::string delta_tau_tau_name = "delta_tau_tau";
const std::string delta_tau_t_name = "delta_tau_t";

// Used by OneSpin QMC

const std::string u_0_t_t_name = "u_0_t_t";
const std::string u_1_t_t_name = "u_1_t_t";
const std::string u_0_tau_tau_name = "u_0_tau_tau";
const std::string u_1_tau_tau_name = "u_1_tau_tau";

const std::string phi_1_t_t_name = "phi_1_t_t";
const std::string phi_1_tau_tau_name = "phi_1_tau_tau";

const std::string delta_up_greater_t_t_name = "delta_up_greater_t_t";
const std::string delta_up_lesser_t_t_name = "delta_up_lesser_t_t";
const std::string delta_up_tau_t_name = "delta_up_tau_t";

const std::string delta_down_greater_t_t_name = "delta_down_greater_t_t";
const std::string delta_down_lesser_t_t_name = "delta_down_lesser_t_t";
const std::string delta_down_tau_t_name = "delta_down_tau_t";

// Used by Determinant QMC

const std::string c_operator_matrices_name = "c_operators_matrices";
const std::string a_operator_matrices_name = "a_operators_matrices";


inline std::string get_filename(const std::string &params_name, const std::string &object_name){
    return params_name + "_" + object_name + "." + extension;
}

