#pragma once

#include "../Measurements.h"

// This measurement is valid only for density-density interactions!

class InitialOccupancies : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        const contour_time *ct_start, *ct_end;
        double tau_end, tau_start;
        double length;
        int n;
        for(int flavor = 0; flavor < conf()->get_n_flavor(); ++flavor){
            // Task: calculate sum of the lengths of antisegments (n = 1) on the imaginary branch
            n = conf()->get_initial_n(flavor); // it as as well the final flavor
            length = 0;
            int k_flavor = conf()->get_flavor_order(flavor);
            DEBUG( std::cout << "Measuring occupancy on Matsubara, flavor " << flavor << std::endl; )
            if(k_flavor == 0){
                // k=0 configuration is either segment or antisegment, occupancy probability = initial_n
                raw_data(flavor) += n * conf()->get_sign().real();
            }
            else{
                for(int i = k_flavor - 1; i >= 0; --i){
                    if(n == 1){
                        // Find the end time of a segment
                        ct_end = &conf()->get_contour_time(creation, flavor, i);
                        // Find the start time of a segment
                        ct_start = &conf()->get_contour_time(annihilation, flavor, i);
                    }
                    else{
                        // Find the end time of an antisegment
                        ct_end = &conf()->get_contour_time(annihilation, flavor, i);
                        // Find the start time of an antisegment
                        ct_start = &conf()->get_contour_time(creation, flavor, i);
                    }
                    if(ct_end->branch != IMAG)
                        break;
                    else
                        tau_end = ct_end->t;
                    if(ct_start->branch != IMAG)
                        tau_start = 0;
                    else
                        tau_start = ct_start->t;
                    assert(tau_end > tau_start);
                    length += (tau_end - tau_start);
                    DEBUG( std::cout << (n == 1 ? "segment" : "antisegment") << " length = " << tau_end - tau_start << std::endl; )
                }
                if(n == 1){
                    // We measured the length of segments, while the length of antisegments is needed
                    length = beta() - length;
                }
                raw_data(flavor) += length / beta() * conf()->get_sign().real();
            }
        }
    }
    void reset() override {
        if(is_master()) {
            mean_data.clear();
            mean_data.resize(conf()->get_n_flavor());
        }
        else{
            raw_data.clear();
            raw_data.resize(conf()->get_n_flavor());
        }
        error_data.clear();
        error_data.resize(conf()->get_n_flavor());
    }
    virtual void save(const std::string &prefix, bool truncate) const {
        const auto suffix = "initial_occupancies";
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + suffix + '_' + std::to_string(rank) + txt_extension,
                          result(), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        const auto suffix = "initial_occupancies";
        result().save(prefix + "_" + suffix + output_extension, output_file_type);
        error().save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        auto res = result();
        auto err = error();
        using std::endl;
        os << name() << ": " << endl;
        for(int flavor = 0; flavor < res.size(); ++flavor) {
            os << " flavor " << flavor << ": n(0) = " << res(flavor) << " +- " << err(flavor) << endl;
        }
    }
    std::string name() const override { return "Local levels initial occupancies"; }
    using data_t = arma::vec; // 1-dim container of real numbers
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
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &global_m) override {
        auto& m = static_cast<InitialOccupancies&>(global_m);
        mpi_all_reduction((N() * result()).eval(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = result() - m.result();
        error_data = N() * error_data % error_data;
        mpi_reduction(error_data, m.error_data);
    }
#endif
};