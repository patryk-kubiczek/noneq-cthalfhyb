#pragma once

#include "../Measurements.h"

class OrderHistogram : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        int k = conf()->get_k();
        if(k >= raw_data.size()){
            raw_data.resize(k + 1);
        }
        ++raw_data(k);
    }
    void reset() override {
        if(is_master()) mean_data.clear();
        else raw_data.clear();
        error_data.clear();
    }
    void save(const std::string &prefix, bool truncate) const override {
        const auto suffix = "order_histogram";
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + suffix + '_' + std::to_string(rank) + txt_extension,
                          result(), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        const auto suffix = "order_histogram";
        result().save(prefix + "_" + suffix + output_extension, output_file_type);
        error().save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        using std::endl;
        os << name() << ": " << endl;
        auto histogram = result();
        auto hist_error = error();
        if(histogram.size() != hist_error.size())
            hist_error.resize(size(histogram));
        for(int k = 0; k < histogram.size(); ++k){
            os << " " << k << ": " << histogram(k) << " +- " << hist_error(k) << endl;
        }
    }
    std::string name() const override { return "Expansion order probability"; }
    using data_t = arma::Col<long long>;
    using result_t = arma::vec;
    result_t result() const {
        return (is_master() ? mean_data : arma::conv_to<result_t>::from(raw_data)) / N();
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
        // Make sure all the histograms have the same size
        auto k_max = raw_data.size();
        MPI_Allreduce(MPI_IN_PLACE, &k_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        raw_data.resize(k_max);
        error_data.resize(k_max);
        auto& m = static_cast<OrderHistogram&>(global_m);
        m.mean_data.resize(k_max);
        m.error_data.resize(k_max);
        mpi_all_reduction((N() * result()).eval(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = result() - m.result();
        error_data = N() * error_data % error_data;
        mpi_reduction(error_data, m.error_data);
    }
#endif
};
