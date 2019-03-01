#pragma once

#include "../local_weight/LocalWeightWrapper.h"
#include "../measurements/Measurements.h"

class QMC {
public:
    QMC(const std::string &name, double beta, double t_max, const InputParams& input_params);
    QMC(const std::string &name, double beta, double t_max, const ModelParams& input_params);
    QMC(const std::string &name, double beta, double t_max,
          const InputParams& input_params, const ModelParams& model_params,
          bool load_delta_up = false, bool load_delta_down = true, bool load_propagator = false);

    // Numerical input
    QMC(const std::string &name, double beta, double t_max,
        double dt, int real_grid_size, double dtau, int imag_grid_size)
            : QMC(name, beta, t_max, InputParams{dt, real_grid_size, dtau, imag_grid_size}) {}

    // Analytical input
    QMC(const std::string &name, double beta, double t_max,
        double eps_0_up, double eps_up, double eps_0_down, double eps_down, double U_0, double U,
        cx_double V_0_up, cx_double V_up,
        double eps_bath_0_up, double eps_bath_up,
        cx_double V_0_down, cx_double V_down,
        double eps_bath_0_down, double eps_bath_down)
            : QMC(name, beta, t_max,
                  ModelParams{eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                              {V_0_up}, {V_up}, {eps_bath_0_up}, {eps_bath_up},
                              {V_0_down}, {V_down}, {eps_bath_0_down}, {eps_bath_down}}) {}

    // Mixed input ("discrete bath")
    QMC(const std::string &name, double beta, double t_max,
        double dt, int real_grid_size, double dtau, int imag_grid_size,
        double eps_0_up, double eps_up, double eps_0_down, double eps_down, double U_0, double U,
        std::vector<cx_double> &V_0_up, std::vector<cx_double> &V_up,
        std::vector<double> &eps_bath_0_up, std::vector<double> &eps_bath_up,
        bool discrete_bath = true)
            : QMC(name, beta, t_max,
                  InputParams{dt, real_grid_size, dtau, imag_grid_size},
                  ModelParams{eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                              V_0_up, V_up, eps_bath_0_up, eps_bath_up,
                              {}, {}, {}, {}}) {}


    void run_measurement(long int N_MC);
    void run_warm_up(long int N_MC);
    void run_measurement(double minutes);
    void run_warm_up(double minutes);

    void reset_measurements() { measurements.reset_measurements(); master.reset_measurements(); }
    void set_n_imag_times(int n) { measurements.set_n_imag_times(n); master.set_n_imag_times(n); }
    void set_n_real_times(int n) { measurements.set_n_real_times(n); master.set_n_real_times(n); }
    void set_n_matsubara(int n) { measurements.set_n_matsubara(n); master.set_n_matsubara(n); }
    void set_block_size(int n) { measurements.set_block_size(n); master.set_block_size(n); }
    void set_n_blocks(int n);
    void print_results();
    void save_results(const std::string &prefix);
    void save_individual_results(const std::string &prefix, bool truncate = false);

#ifdef USE_MPI
    void collect_results() { reduce_measurements(); }
#endif

    friend std::ostream &operator<<(std::ostream &os, const QMC &qmc);

private:

#ifdef USE_MPI
    #ifndef USE_MPI4PY
    MPI_Initializer mpi_initializer;
    #endif
    void reduce_measurements();
    void reduce_exception_counter();
#endif

    local_weight_type lw_type;
    bath_weight_type bw_type;

    RandomNumberGenerator random_number_generator;
    Model model;
    Configuration configuration;
    LocalWeightWrapper<local_weight_t> local_weight;
    BathWeight bath_weight;
    Moves moves;
    Measurements measurements;
    Measurements master;

    bool is_measuring = false;

    void initialize_configuration_and_weights();
    void initialize_random_number_generator();
    void perform_move(int flavor);
    void run(long int N_MC);
    void run(double minutes);

    void print_state(std::ostream &os) const;
    void crosscheck_move();

};