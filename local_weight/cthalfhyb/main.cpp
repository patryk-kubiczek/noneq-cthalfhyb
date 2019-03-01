#include "../../model/cthalfhyb/Model.h"
#include "LocalWeight.h"
#include "../../discrete_bath/DiscreteBathWeight.h"

int main(){
    std::string params_name = "../../../input/test";

    double beta = 1;
    double t_max = 1;

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

    int real_grid_size = 250;
    int imag_grid_size = 250;

    InputParams input_params = {1. / (real_grid_size - 1), real_grid_size, beta / (imag_grid_size - 1) , imag_grid_size};

    ModelParams model_params = {eps_0_up, eps_up, eps_0_down, eps_up, U_0, U,
                                V_0_up, V_up, eps_bath_0_up, eps_bath_up, V_0_down, V_down, eps_bath_0_down, eps_bath_down};

    Model model(params_name, beta, t_max, input_params, model_params, true, true, true);
    std::cout << model << std::endl;

    RandomNumberGenerator rng;
    int seed = std::random_device{}();
    //std::cin >> seed;
    rng.initialize(seed);

    Configuration conf(model.t_max, model.beta, model.n_flavor);
    conf.randomly_generate(1, rng);

    std::cout << conf << std::endl;

    LocalWeight local_weight(model, conf, NUMERIC_LW);

    std::cout << local_weight << std::endl;

    LocalWeight analytic_local_weight(model, conf, ANALYTIC_LW);

    std::cout << analytic_local_weight << std::endl;

    DiscreteBathWeight discrete_bath(model, conf, DISCRETE_BATH);

    std::cout << discrete_bath << std::endl;

    int n_1, n_2, n_3;

    std::cout << local_weight.density_matrix_at_t_max(&n_1) << std::endl;
    std::cout << analytic_local_weight.density_matrix_at_t_max(&n_2) << std::endl;
    std::cout << discrete_bath.density_matrix_at_t_max(&n_3) << std::endl;

    std::cout << "n_down at t_max: " << n_1 << ", " << n_2 << ", " << n_3 << std::endl;
}