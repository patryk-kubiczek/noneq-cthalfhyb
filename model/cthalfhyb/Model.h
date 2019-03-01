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
    Model(const std::string &name, double beta, double t_max, const InputParams& input_params);
    Model(const std::string &name, double beta, double t_max, const ModelParams& model_params);
    Model(const std::string &name, double beta, double t_max,
          const InputParams& input_params, const ModelParams& model_params,
          bool load_delta_up, bool load_delta_down, bool load_propagator);

    const std::string params_name;
    const double beta;
    const double t_max;

    InputParams input_params;
    ModelParams model_params;
    int n_bath_up;

    const int n_flavor = 1;

    const ContourPropagator u_0;
    const ContourPropagator u_1;
    const ContourPhi phi_1;
    const ContourDeltaSingle delta_up;
    const ContourDeltaSingle delta_down;

    const AnalyticU<0> u_0_analytic;
    const AnalyticU<1> u_1_analytic;
    const AnalyticPhi1 phi_1_analytic;
    const AnalyticDelta delta_up_analytic;
    const AnalyticDelta delta_down_analytic;

    friend std::ostream &operator<<(std::ostream &os, const Model &model);

private:
    Model(const std::string &name, double beta, double t_max);

    cx_tetracube u_0_t_t_data;
    cx_tetracube u_1_t_t_data;
    arma::cube u_0_tau_tau_data;
    arma::cube u_1_tau_tau_data;

    arma::cx_mat phi_1_t_t_data;
    arma::vec phi_1_tau_tau_data;

    arma::cx_mat delta_up_greater_t_t_data;
    arma::cx_mat delta_up_lesser_t_t_data;
    arma::cx_mat delta_up_tau_t_data;
    arma::cx_mat delta_down_greater_t_t_data;
    arma::cx_mat delta_down_lesser_t_t_data;
    arma::cx_mat delta_down_tau_t_data;

    void load_delta_up_data();
    void load_delta_down_data();
    void load_propagator_data();

    bool loaded_delta_up = false;
    bool loaded_delta_down = false;
    bool loaded_propagator = false;
};