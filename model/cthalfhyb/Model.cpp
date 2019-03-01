#include "Model.h"

Model::Model(const std::string &name, double beta, double t_max)
        : params_name(name), beta(beta), t_max(t_max),
          u_0(u_0_t_t_data, u_0_tau_tau_data, input_params.dt, input_params.dtau, beta),
          u_1(u_1_t_t_data, u_1_tau_tau_data, input_params.dt, input_params.dtau, beta),
          phi_1(phi_1_t_t_data, phi_1_tau_tau_data, input_params.dt, input_params.dtau, beta),
          delta_up(delta_up_greater_t_t_data, delta_up_lesser_t_t_data, delta_up_tau_t_data, input_params.dt, input_params.dtau, beta),
          delta_down(delta_down_greater_t_t_data, delta_down_lesser_t_t_data, delta_down_tau_t_data, input_params.dt, input_params.dtau, beta),
          u_0_analytic(beta, model_params),
          u_1_analytic(beta, model_params),
          phi_1_analytic(beta, model_params),
          delta_up_analytic(beta, model_params, 0),
          delta_down_analytic(beta, model_params, 1) {}

Model::Model(const std::string &name, double beta, double t_max, const InputParams &params)
        : Model(name, beta, t_max) {
    input_params = params;
    load_delta_down_data();
    load_propagator_data();
    n_bath_up = u_0_tau_tau_data.n_rows - 1;
}

Model::Model(const std::string &name, double beta, double t_max, const ModelParams &params)
        : Model(name, beta, t_max) {
    model_params = params;
    n_bath_up = model_params.eps_bath_0_up.size();
}

Model::Model(const std::string &name, double beta, double t_max, const InputParams &input_params,
             const ModelParams &model_params, bool load_delta_up, bool load_delta_down, bool load_propagator)
        : Model(name, beta, t_max)  {
    this->input_params = input_params;
    this->model_params = model_params;
    if(load_propagator) {
        load_propagator_data();
        n_bath_up = u_0_tau_tau_data.n_rows - 1;
    }
    else
        n_bath_up = model_params.eps_bath_0_up.size();
    if(load_delta_up) load_delta_up_data();
    if(load_delta_down) load_delta_down_data();
}

void Model::load_delta_up_data() {
    delta_up_greater_t_t_data.load(get_filename(params_name, delta_up_greater_t_t_name));
    delta_up_lesser_t_t_data.load(get_filename(params_name, delta_up_lesser_t_t_name));
    delta_up_tau_t_data.load(get_filename(params_name, delta_up_tau_t_name));
    loaded_delta_up = true;
}

void Model::load_delta_down_data() {
    delta_down_greater_t_t_data.load(get_filename(params_name, delta_down_greater_t_t_name));
    delta_down_lesser_t_t_data.load(get_filename(params_name, delta_down_lesser_t_t_name));
    delta_down_tau_t_data.load(get_filename(params_name, delta_down_tau_t_name));
    loaded_delta_down = true;
}

void Model::load_propagator_data() {
    u_0_t_t_data.load(get_filename(params_name, u_0_t_t_name), input_params.real_grid_size, input_params.real_grid_size);
    u_1_t_t_data.load(get_filename(params_name, u_1_t_t_name), input_params.real_grid_size, input_params.real_grid_size);
    u_0_tau_tau_data.load(get_filename(params_name, u_0_tau_tau_name));
    u_1_tau_tau_data.load(get_filename(params_name, u_1_tau_tau_name));
    phi_1_t_t_data.load(get_filename(params_name, phi_1_t_t_name));
    phi_1_tau_tau_data.load(get_filename(params_name, phi_1_tau_tau_name));
    loaded_propagator = true;
}


std::ostream &operator<<(std::ostream &os, const Model &model) {
    using std::endl;
    os << "*** Model data ***" << endl;
    os << "Model name: " << model.params_name << endl;
    os << "beta: " << model.beta << endl;
    os << "t_max: " << model.t_max << endl;
    os << "n_bath_up: " << model.n_bath_up << endl;
    if(model.loaded_delta_down || model.loaded_delta_up || model.loaded_propagator) {
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

