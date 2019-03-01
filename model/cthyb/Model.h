#pragma once

#include <ostream>
#include "../../auxiliary_functions/constants.h"
#include "../../data_types/Tetracube.h"
#include "../ContourDelta.h"
#include "../ContourPropagator.h"
#include "../analytic_expressions.h"
#include "../params.h"

class Model{
public:
    Model(const std::string &name, double beta, double t_max);
    Model(const std::string &name, double beta, double t_max, const InputParams& input_params);
    Model(const std::string &name, double beta, double t_max, const ModelParams& model_params);

    const std::string params_name;
    const double beta;
    const double t_max;

    InputParams input_params;
    ModelParams model_params;

    int n_flavor;
    const int n_bath_up = 0;

    const ContourPropagator bare_propagator;
    const ContourDelta delta;
    
    const AnalyticP bare_propagator_analytic;
    const AnalyticDelta delta_up_analytic;
    const AnalyticDelta delta_down_analytic;

    arma::cx_cube c_operator_matrices;
    arma::cx_cube a_operator_matrices;

    friend std::ostream &operator<<(std::ostream &os, const Model &model);

private:
    cx_tetracube p_t_t_data;
    arma::cube p_tau_tau_data;

    arma::field<arma::cx_mat> delta_greater_t_t_data;
    arma::field<arma::cx_mat> delta_lesser_t_t_data;
    arma::field<arma::cx_mat> delta_tau_t_data;
    
    void load_delta_data();
    void load_propagator_data();
    void load_operator_matrices();

    bool loaded_delta = false;
    bool loaded_propagator = false;
};