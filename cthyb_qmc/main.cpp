#include "QMC.h"

int main() {
    std::string output_name = "output/CTHYB/test";
    std::string params_name = "../../input/test";


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

    ModelParams model_params = {eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                                 V_0_up, V_up, eps_bath_0_up, eps_bath_up, V_0_down, V_down, eps_bath_0_down, eps_bath_down};

    double input_t_max = 1;
    int real_grid_size = 250;
    int imag_grid_size = 250;
    double dt = input_t_max / (real_grid_size - 1);
    double dtau = beta / (imag_grid_size - 1);

    InputParams input_params = {dt, real_grid_size, dtau, imag_grid_size};

    //QMC qmc(params_name, beta, t_max, model_params);
    QMC qmc(params_name, beta, t_max, input_params);


    long int N_warm_up = 1000;
    long int N_MC =      100000;
    double time_MC = 1;
    double time_warm_up = 0.05 * time_MC;

    int n_real = 10;
    int n_imag = 50;
    int n_block = 1;

    qmc.set_n_real_times(n_real);
    qmc.set_n_imag_times(n_imag);
    qmc.set_block_size(n_block);

    qmc.reset_measurements();
    qmc.run_warm_up(N_warm_up);
    qmc.reset_measurements();
    qmc.run_measurement(N_MC);
#ifdef USE_MPI
    qmc.collect_results();
#endif
    qmc.print_results();
    qmc.save_results(output_name);


}
