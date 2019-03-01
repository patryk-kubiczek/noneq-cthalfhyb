#include "Model.h"

Model::Model(const std::string &name, double beta, double t_max)
        : params_name(name), beta(beta), t_max(t_max),
          bare_propagator(p_t_t_data, p_tau_tau_data, input_params.dt, input_params.dtau, beta),
          delta(delta_greater_t_t_data, delta_lesser_t_t_data, delta_tau_t_data, input_params.dt, input_params.dtau, beta),
          bare_propagator_analytic(beta, model_params),
          delta_up_analytic(beta, model_params, 0),
          delta_down_analytic(beta, model_params, 1) {
    load_operator_matrices();
}

Model::Model(const std::string &name, double beta, double t_max, const InputParams &params)
        : Model(name, beta, t_max) {
    input_params = params;
    load_delta_data();
    load_propagator_data();
    n_flavor = delta_greater_t_t_data.n_rows;
}

Model::Model(const std::string &name, double beta, double t_max, const ModelParams &params)
        : Model(name, beta, t_max) {
    model_params = params;
    n_flavor = 2;
}

void Model::load_delta_data() {
    delta_greater_t_t_data.load(get_filename(params_name, delta_greater_t_t_name));
    delta_lesser_t_t_data.load(get_filename(params_name, delta_lesser_t_t_name));
    delta_tau_t_data.load(get_filename(params_name, delta_tau_t_name));
    loaded_delta = true;
}



void Model::load_propagator_data() {
    p_t_t_data.load(get_filename(params_name, p_t_t_name), input_params.real_grid_size, input_params.real_grid_size);
    p_tau_tau_data.load(get_filename(params_name, p_tau_tau_name));
    loaded_propagator = true;
}

void Model::load_operator_matrices() {
    c_operator_matrices.load(get_filename(params_name, c_operator_matrices_name));
    a_operator_matrices.load(get_filename(params_name, a_operator_matrices_name));
}

std::ostream &operator<<(std::ostream &os, const Model &model) {
    using std::endl;
    os << "*** Model data ***" << endl;
    os << "Model name: " << model.params_name << endl;
    os << "beta: " << model.beta << endl;
    os << "t_max: " << model.t_max << endl;
    if(model.loaded_delta || model.loaded_propagator) {
        os << "Numerical input parameters:" << endl;
        os << "dtau: " << model.input_params.dtau << ", dt: " << model.input_params.dt << endl;
    }
    if(!model.loaded_propagator) {
        os << "Model input parameters:" << endl;
        os << "eps_0_up: " << model.model_params.eps_0_up << ", eps_up: " << model.model_params.eps_up << endl;
        os << "eps_0_down: " << model.model_params.eps_0_down << ", eps_down: " << model.model_params.eps_down << endl;
        os << "U_0: " << model.model_params.U_0 << ", U: " << model.model_params.U << endl;
    }
    return os;
}
