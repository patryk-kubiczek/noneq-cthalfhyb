#pragma once

#include <unordered_map>
#include <ostream>
#include <cstdio>
#include <type_traits>

#include "../auxiliary_functions/save_functions.h"
#include "../auxiliary_functions/element_wise_functions.h"
#include "../auxiliary_functions/complex_numbers.h"

#include "../configuration/Configuration.h"
#include "../local_weight/LocalWeightWrapper.h"
#include "../moves/Moves.h"

#ifdef USE_MPI
#include "../auxiliary_functions/mpi_wrappers.h"
#endif



const auto output_file_type = arma::raw_ascii;
const std::string output_extension =  ".out";
const std::string txt_extension =  ".txt";


class Measurements {
public:
    Measurements(const Model &model,
                 const Configuration &conf,
                 const local_weight_t &local_weight,
                 const BathWeight &bath_weight,
                 const Moves &moves,
                 bool is_master = false)
            : t_max(model.t_max), beta(model.beta),
              n_bath(model.n_bath_up),
              conf(&conf), local_weight(&local_weight),
              bath_weight(&bath_weight), moves(&moves),
              is_master(is_master){
        initialize_measurements();
    }

    class MeasurementBase{
    public:
        explicit MeasurementBase(const Measurements &measurements) : meas(&measurements) {}
        virtual void measure() = 0;
        virtual void reset() = 0;
        virtual void save(const std::string &prefix, bool truncate) const = 0;
        virtual void master_save(const std::string &prefix) const = 0;
        virtual void print_normalized(std::ostream& os) const = 0;
        virtual std::string name() const = 0;
        virtual ~MeasurementBase() = default;
    protected:
        const Configuration* conf() const { return meas->conf; }
        const local_weight_t* local_weight() const { return meas->local_weight; }
        const BathWeight* bath_weight() const { return meas->bath_weight; }
        long long N() const { return meas->measurement_counter; }
        int n_times() const { return meas->n_times; }
        int n_imag_times() const { return meas->n_imag_times; }
        int n_matsubara() const { return meas->n_matsubara; }
        double t_max() const { return meas->t_max; }
        double beta() const { return meas->beta; }
        int n_bath() const { return meas->n_bath; }
        bool is_master() const { return meas->is_master; }
        cx_double mean_sign() const { return meas->mean_sign(); }
    private:
        const Measurements *meas;
    public:
#ifdef USE_MPI
        virtual void reduce(MeasurementBase &master_measurement) = 0;
#endif
    };

    void increment_step_counter() { ++step_counter; }
    auto get_step_counter() { return step_counter; }

    void increment_accepted_moves_counter(int move_number) { ++acceptance_rates[move_number]; }

    void increment_exception_counter(const std::string &exception_name);

    void do_measurements();

    void reset_measurements();

    void set_n_imag_times(int n) { n_imag_times = n; }
    void set_n_real_times(int n) { n_times = n; }
    void set_n_matsubara(int n) { n_matsubara = n; }
    void set_block_size(int n) { block_size = n; }


    friend std::ostream &operator<<(std::ostream &os, const Measurements &meas);

    MeasurementBase& measurement(int i) { return *(measurements[i]); }
    MeasurementBase& measurement(int i) const { return *(measurements[i]); }
    auto measurements_size() const { return measurements.size(); }


    auto& step_counter_data() { return step_counter; }
    auto& measurement_counter_data() { return measurement_counter; }
    auto& exception_counter_data() { return exception_counter; }
    auto& acceptance_rates_data() { return acceptance_rates; }

private:
    double t_max;
    double beta;
    int n_bath;

    const Configuration *conf;
    const local_weight_t *local_weight;
    const BathWeight *bath_weight;
    const Moves *moves;

    long long step_counter = 0;
    long long measurement_counter = 0;
    std::unordered_map<std::string, int> exception_counter;
    std::vector<long long> acceptance_rates;

    std::vector<std::unique_ptr<MeasurementBase>> measurements;

    int n_times = 10;
    int n_imag_times = 50;
    int n_matsubara = 50;
    int block_size = 10;

    bool is_master = false;

    void initialize_measurements();

    cx_double mean_sign() const;
    double mean_acceptance() const;
    void print_acceptance_rates(std::ostream &os) const;
};