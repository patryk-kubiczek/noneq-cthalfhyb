#include <chrono>

#include "QMC.h"


QMC::QMC(const std::string &name, double beta, double t_max, const InputParams &input_params)
        : bw_type(NUMERIC_BW), lw_type(NUMERIC_LW),
          model(name, beta, t_max, input_params),
          configuration(t_max, beta, model.n_flavor),
          bath_weight(model, configuration, NUMERIC_BW),
          local_weight(model, configuration, NUMERIC_LW),
          moves(model, configuration, local_weight.get(), bath_weight, random_number_generator),
          measurements(model, configuration, local_weight.get(), bath_weight, moves),
          master(model, configuration, local_weight.get(), bath_weight, moves, true) {

    initialize_random_number_generator();
    initialize_configuration_and_weights();
    std::cout << (*this) << std::endl;
}

QMC::QMC(const std::string &name, double beta, double t_max, const ModelParams &model_params)
        : bw_type(ANALYTIC_BW), lw_type(ANALYTIC_LW),
          model(name, beta, t_max, model_params),
          configuration(t_max, beta, model.n_flavor),
          bath_weight(model, configuration, ANALYTIC_BW),
          local_weight(model, configuration, ANALYTIC_LW),
          moves(model, configuration, local_weight.get(), bath_weight, random_number_generator),
          measurements(model, configuration, local_weight.get(), bath_weight, moves),
          master(model, configuration, local_weight.get(), bath_weight, moves, true) {
    initialize_random_number_generator();
    initialize_configuration_and_weights();
    std::cout << (*this) << std::endl;
}


void QMC::initialize_configuration_and_weights() {
    // Generating the initial configuration. Default: k = 2 configuration
    configuration.randomly_generate(2, random_number_generator);
    local_weight.get().generate_from_conf();
    bath_weight.generate_from_conf();
    configuration.set_weight(local_weight.get().get_w_loc() * cx_long_double{bath_weight.get_w_bath()});
}

void QMC::run_measurement(long int N_MC) {
#ifdef USE_MPI
    if(mpi_rank() == 0)
#endif
        std::cout << std::endl << "STARTING MEASUREMENT" << std::endl << std::endl;
    if(!is_measuring){
        //measurements.reset_measurements();
        is_measuring = true;
    }
    run(N_MC);
}

void QMC::run_warm_up(long int N_MC) {
#ifdef USE_MPI
    if(mpi_rank() == 0)
#endif
        std::cout << std::endl << "STARTING WARM UP" << std::endl << std::endl;
    if(is_measuring){
        //measurements.reset_measurements();
        is_measuring = false;
    }
    run(N_MC);
}

void QMC::run_measurement(double minutes) {
#ifdef USE_MPI
    if(mpi_rank() == 0)
#endif
        std::cout << std::endl << "STARTING MEASUREMENT" << std::endl << std::endl;
    if(!is_measuring){
        //measurements.reset_measurements();
        is_measuring = true;
    }
    run(minutes);
}

void QMC::run_warm_up(double minutes) {
#ifdef USE_MPI
    if(mpi_rank() == 0)
#endif
        std::cout << std::endl << "STARTING WARM UP" << std::endl << std::endl;
    if(is_measuring){
        //measurements.reset_measurements();
        is_measuring = false;
    }
    run(minutes);
}


void QMC::print_results() {
#ifndef USE_MPI
    std::cout << measurements << std::endl;
#endif
#ifdef USE_MPI
    if(mpi_rank() == 0){
        std::cout << master << std::endl;
    }
#endif
}

void QMC::save_results(const std::string &prefix) {
#ifndef USE_MPI
    const auto &m = measurements;
#endif
#ifdef USE_MPI
    const auto &m = master;
    if(mpi_rank() == 0){
#endif
    std::ofstream output_file(prefix + "_qmc.out");
    if(output_file){
        output_file << model << std::endl;
        output_file << m << std::endl;
        output_file.close();
    }
    for(int i = 0; i < m.measurements_size(); ++i)
        m.measurement(i).master_save(prefix);
#ifdef USE_MPI
    }
#endif
}

void QMC::save_individual_results(const std::string &prefix, bool truncate) {
#ifndef USE_MPI
    const auto &m = measurements;
#endif
#ifdef USE_MPI
    const auto &m = master;
    if(mpi_rank() == 0){
#endif
    std::ofstream output_file(prefix + "_qmc.txt");
    if(output_file){
        output_file << model << std::endl;
        output_file << m << std::endl;
        output_file.close();
    }
#ifdef USE_MPI
    }
#endif
    for(int i = 0; i < measurements.measurements_size(); ++i)
        measurements.measurement(i).save(prefix, truncate);

}

void QMC::perform_move(int flavor) {
    bool is_move_accepted;
    auto &move = moves.random_move();
    try{
        DEBUG( std::cout << "STEP no. " << measurements.get_step_counter() << ":" << std::endl;
                       std::cout << "Move [" << move.name() << ", flavor: " << flavor << "] " << std::endl; )
        is_move_accepted = move.execute(flavor);
        DEBUG(
        if(is_move_accepted){
            std::cout << "ACCEPTED" << std::endl;
            crosscheck_move();
            print_state(std::cout);
            std::cout << std::endl;
        }
        else
            std::cout << "REJECTED" << std::endl << std::endl;
            print_state(std::cout);
            std::cout << std::endl;
        )
        if(is_measuring) {
            measurements.do_measurements();
            if(is_move_accepted)
                measurements.increment_accepted_moves_counter(move.index);
        }
        measurements.increment_step_counter();
    }
    catch (const std::runtime_error &error){
        measurements.increment_exception_counter(error.what());
    }
}

void QMC::run(long int N_MC) {
    int flavor = 0;
    for(long int i = 0; i < N_MC; ++i){
        flavor = random_number_generator.random_int(0, model.n_flavor - 1);
        perform_move(flavor);
    }
}

void QMC::run(double minutes) {
    using namespace std::chrono;
    const auto clock_begin = steady_clock::now();
    auto duration_in_min = [&clock_begin](){
        return duration<double, std::ratio<60>>(steady_clock::now() - clock_begin).count(); };

    int flavor = 0;
    const int n_moves_before_time_check = 100;
    while(duration_in_min() < minutes){
        for(int i = 0; i < n_moves_before_time_check; ++i){
            flavor = random_number_generator.random_int(0, model.n_flavor - 1);
            perform_move(flavor);
        }
    }
}

void QMC::print_state(std::ostream &os) const {
    using std::endl;
    os << configuration << endl;
    os << local_weight.get() << endl;
    os << bath_weight << endl;
}

void QMC::crosscheck_move() {
    const double epsilon = 0.00001;

    Configuration new_conf(configuration);

    LocalWeightWrapper<local_weight_t> new_local_weight(model, new_conf, lw_type);
    new_local_weight.get().generate_from_conf();

    BathWeight new_bath_weight(model, new_conf, bw_type);
    new_bath_weight.generate_from_conf();

    new_conf.set_weight(new_local_weight.get().get_w_loc() * cx_long_double{new_bath_weight.get_w_bath()});

    auto local_weight_from_move = local_weight.get().get_w_loc();
    auto local_weight_from_scratch = new_local_weight.get().get_w_loc();
    std::cout << "w_loc from move: " << local_weight_from_move << ", w_loc from scratch: " << local_weight_from_scratch << std::endl;

    auto bath_weight_from_move = bath_weight.get_w_bath();
    auto bath_weight_from_scratch = new_bath_weight.get_w_bath();
    std::cout << "w_b from move: " << bath_weight_from_move << ", w_b from scratch: " << bath_weight_from_scratch << std::endl;

    auto abs_w_from_move = configuration.get_abs_weight();
    auto abs_w_from_scratch = new_conf.get_abs_weight();
    //std::cout << "matrix element from move: " << local_weight.get_contributing_matrix_element()
    //          << ", matrix element from scratch: " << new_local_weight.get_contributing_matrix_element() << std::endl;
    std::cout << "|w_C| from move: " << abs_w_from_move << ", |w_C| from scratch: " << abs_w_from_scratch << std::endl;

    auto phase_w_from_move = std::arg(configuration.get_sign());
    auto phase_w_from_scratch = std::arg(new_conf.get_sign());
    std::cout << "phase from move: " << phase_w_from_move << ", phase from scratch: " << phase_w_from_scratch << std::endl;

    if(std::abs((abs_w_from_move - abs_w_from_scratch) / abs_w_from_move) > epsilon){
        std::cout << "|w_C| not updated correctly" << std::endl;
        std::cout << configuration;
        std::cout << new_conf;
        //assert(false);
        throw std::runtime_error("|w_C| not updated correctly, diff: " + std::to_string(abs_w_from_move - abs_w_from_scratch));
    }
    if(std::abs(phase_w_from_move - phase_w_from_scratch) > epsilon){
        if(std::abs(std::abs(phase_w_from_move - phase_w_from_scratch) - 2 * PI) > epsilon){
            std::cout << "phase not updated correctly" << std::endl;
            throw std::runtime_error("phase not updated correctly, diff: " + std::to_string(phase_w_from_move - phase_w_from_scratch));
        }
    }
}

void QMC::initialize_random_number_generator() {
    random_number_generator.initialize(std::random_device{}());
    #ifdef USE_MPI
    unsigned int seed;
    if(mpi_rank() == 0)
        seed = std::random_device{}();
    mpi_broadcast(seed);
    random_number_generator.initialize(seed, mpi_size(), mpi_rank());
    #endif
}


std::ostream &operator<<(std::ostream &os, const QMC &qmc) {
#ifdef USE_MPI
    if(mpi_rank() == 0){
#endif
    os << "*** Non-equilibrium CT-HYB QMC ***" << std::endl;
    os << std::endl;
    os << qmc.model;
    os << qmc.moves << std::endl;
#ifdef USE_MPI
    os << "Number of parallel processes: " << mpi_size() << std::endl << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    DEBUG( std::cout << "Current state:" << std::endl;
                   qmc.print_state(os); )
    return os;
}




#ifdef USE_MPI
void QMC::reduce_measurements() {

    auto& local_m = measurements;
    auto& global_m = master;

    mpi_reduction(local_m.step_counter_data(), global_m.step_counter_data());
    mpi_all_reduction(local_m.measurement_counter_data(), global_m.measurement_counter_data());
    mpi_reduction(local_m.acceptance_rates_data(), global_m.acceptance_rates_data());
    reduce_exception_counter();

    for (int i = 0; i < local_m.measurements_size(); ++i){
        local_m.measurement(i).reduce(global_m.measurement(i));
    }
}

void QMC::reduce_exception_counter() {

    auto& local_m = measurements;
    auto& global_m = master;

    int n_processes = mpi_size();
    for(int i = 0; i < n_processes; ++i){
        for(const auto& exception : local_m.exception_counter_data()){
            if(mpi_rank() == i){
                // Send exception name and exception counter from i to 0
                MPI_Send(exception.first.c_str(), exception.first.length() + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
                MPI_Send(&(exception.second), 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            }
            if(mpi_rank() == 0){
                // Receive exception name at 0 from i as char then convert it to string
                MPI_Status status;
                MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
                int char_length;
                MPI_Get_count(&status, MPI_CHAR, &char_length);
                char *received_char = new char[char_length];
                MPI_Recv(received_char, char_length, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
                std::string exception_name(received_char);
                delete[](received_char);
                // Receive exception counter at 0 from i
                int exception_count;
                MPI_Recv(&exception_count, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                // Update master exception counter
                if(global_m.exception_counter_data().count(exception_name) == 0)
                    global_m.exception_counter_data()[exception_name] = exception_count;
                else
                    global_m.exception_counter_data()[exception_name] += exception_count;
            }
        }
    }
}
#endif