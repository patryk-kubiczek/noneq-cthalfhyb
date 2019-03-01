#include "Model.h"

int main(){
    std::string params_name = "../../../input/test";

    double beta = 1;
    double t_max = 0.5;

    double eps_0_up = 1;
    double eps_up = -1;
    double eps_0_down = 1;
    double eps_down = -1;

    double U_0 = 3;
    double U = 3;

    std::vector<cx_double> V_0_up = {1};
    std::vector<cx_double> V_up = {1};
    std::vector<double> eps_bath_0_up = {0};
    std::vector<double> eps_bath_up = {0};
    std::vector<cx_double> V_0_down = {1};
    std::vector<cx_double> V_down = {1};
    std::vector<double> eps_bath_0_down = {0};
    std::vector<double> eps_bath_down = {0};

    ModelParams model_params =  {eps_0_up, eps_up, eps_0_down, eps_0_up, U_0, U,
                                 V_0_up, V_up, eps_bath_0_up, eps_bath_up, V_0_down, V_down, eps_bath_0_down, eps_bath_down};

    double input_t_max = 1;
    int real_grid_size = 100;
    int imag_grid_size = 100;
    double dt = input_t_max / (real_grid_size - 1);
    double dtau = beta / (imag_grid_size - 1);

    InputParams input_params = {dt, real_grid_size, dtau, imag_grid_size};

    Model model_full_numeric(params_name, beta, t_max, input_params);
    std::cout << model_full_numeric;
    //std::cout << model_full_numeric.u_1({IMAG, 0.5}, {IMAG, 0.2});

    Model model_full_analytic (params_name, beta, t_max, model_params);
    std::cout << model_full_analytic;
    
}

