#pragma once

#include "../Measurements.h"

class MeanOrder : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;
    void measure() override {
        raw_data += conf()->get_k();
    }
    virtual void reset() override {
        if(is_master()) mean_data = 0;
        else raw_data = 0;
        error_data = 0;
    }
    virtual void save(const std::string &prefix, bool truncate) const override {
        const auto suffix = "mean_order";
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + suffix + '_' + std::to_string(rank) + txt_extension,
                          result(), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        arma::Col<result_t> res = {result()};
        arma::Col<result_t> err = {error()};
        const auto suffix = "mean_order";
        res.save(prefix + "_" + suffix + output_extension, output_file_type);
        err.save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        os << name() << ": " << result() << " +- " << error() << std::endl;
    }
    std::string name() const override { return "Mean expansion order"; }
    using data_t = long long;
    using result_t = double;
    result_t result() const { return (is_master() ? mean_data : result_t(raw_data)) / N(); }
    result_t error() const { return std::sqrt(error_data / N()); }
private:
    union {
        data_t raw_data;
        result_t mean_data;
    };
    result_t error_data;
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &master_m) override {
        auto& m = static_cast<MeanOrder&>(master_m);
        mpi_all_reduction(N() * result(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = (result() - m.result());
        error_data *= N() * error_data;
        mpi_reduction(error_data, m.error_data);
    }
#endif
};
