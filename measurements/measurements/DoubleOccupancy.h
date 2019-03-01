#pragma once

#include "../Measurements.h"

// This measurement is valid only for density-density interactions and Determinant QMC!
// Provides WRONG results for t_max = 0 - equilibrium calculation not yet implemented

class DoubleOccupancy : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    arma::rowvec occupancies(int flavor) {
        arma::rowvec result(2 * n_times(), arma::fill::zeros);
        const int k = conf()->get_flavor_order(flavor);

        if(k == 0){
            for(int i = 0; i < 2 * n_times(); ++i){
                result(i) = conf()->get_initial_n(flavor);
            }
        }
        else{
            auto segments = conf()->get_segments(flavor);
            int j = 2 * k - 1;
            // Here we are traversing the PLUS branch of the contour t = 0 -> t = t_max
            for(int i = 0; i < n_times(); ++i){
                while (!segments[j].contains(contour_time{PLUS, i * dt})) j = (j + 1) % (2 * k);
                result(i) += segments[j].n();
            }
            // Here we are traversing the MINUS branch of the contour t = t_max -> t = 0
            for(int i = n_times() - 1; i >= 0; --i){
                while (!segments[j].contains(contour_time{MINUS, i * dt})) j = (j + 1) % (2 * k);
                result(n_times() + i) += segments[j].n();
            }
        }
        return result;
    }

    void measure() override {
        arma::rowvec double_occupancy = occupancies(0) % occupancies(1);
        raw_data += (double_occupancy.head(n_times()) + double_occupancy.tail(n_times())) / 2. * conf()->get_sign().real();

    }
    void reset() override {
        dt = t_max() / (n_times() - 1);
        if(is_master()) {
            mean_data.zeros(1, n_times());
        }
        else{
            raw_data.zeros(1, n_times());
        }
        error_data.zeros(1, n_times());
    }
    void save(const std::string &prefix, bool truncate) const override {
        const auto suffix = "double_occupancy";
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + suffix + '_' + std::to_string(rank) + txt_extension,
                          result(), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        const auto suffix = "double_occupancy";
        result().save(prefix + "_" + suffix + output_extension, output_file_type);
        error().save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        auto res = result();
        auto err = error();
        using std::endl;
        os << name() << ": " << endl;
        os << " d(t) = [";
        bool first = true;
        for(int i : std::vector<int>{0, n_times() / 4, n_times() / 2, 3 * n_times() / 4, n_times() - 1}){
            if(first) first = false;
            else os << ", ";
            os << res(i) << " +- " << err(i);
        }
        os << "]" << endl;

    }
    std::string name() const override { return "Time-dependent double occupancy"; }
    using data_t = arma::mat;
    using result_t = data_t;
    result_t result() const {
        return (is_master() ? mean_data : raw_data / mean_sign().real()) /  N();
    }
    result_t error() const {
        return arma::sqrt(error_data / N());
    }
private:
    data_t raw_data;
    result_t mean_data;
    result_t error_data;
    double dt;
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &global_m) override {
        auto& m = static_cast<DoubleOccupancy&>(global_m);
        mpi_all_reduction((N() * result()).eval(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = result() - m.result();
        error_data = N() * error_data % error_data;
        mpi_reduction(error_data, m.error_data);
    }
#endif
};
