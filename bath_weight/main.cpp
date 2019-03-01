#include "BathWeight.h"

int main(){
    std::string params_name = "../../../input/test";

    int n_flavor = 2;
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

    Model model(params_name, beta, t_max, model_params);
    std::cout << model << std::endl;

    RandomNumberGenerator rng;

    Configuration conf(t_max, beta, n_flavor);
    conf.randomly_generate(3, rng);

    std::cout << conf;

    BathWeight bath_weight(model, conf, ANALYTIC_BW);

    bath_weight.generate_from_conf();

    std::cout << bath_weight.get_w_bath() << std::endl;

    std::cout << bath_weight;
}

