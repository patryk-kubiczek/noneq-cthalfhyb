#pragma once

#include "../Measurements.h"

class MeanComplexSign : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        assert(std::abs(conf()->get_sign()) - 1.0 < 0.00001);
        if(!conf()->has_operators_on_real_branch())
            raw_data += conf()->get_sign();
    }

    void reset() override {
        if(is_master()) mean_data = 0.;
        else raw_data = 0.;
        error_data = 0.;
    }

    void save(const std::string &prefix, bool truncate) const override {
        const auto suffix = "mean_complex_sign";
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
        const auto suffix = "mean_complex_sign";
        res.save(prefix + "_" + suffix + output_extension, output_file_type);
        err.save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }

    void print_normalized(std::ostream &os) const override {
        os << name() << ": ";
        print_complex_number(os, result(), error());
        os << std::endl;
    }

    std::string name() const override { return "Mean complex sign"; }

    using data_t = cx_double;
    using result_t = data_t;
    result_t result() const { return (is_master() ? mean_data : result_t(raw_data)) / result_t::value_type(N()); }
    result_t error() const { return element_wise::sqrt(error_data / result_t::value_type(N())); }
private:
    data_t raw_data;
    result_t mean_data;
    result_t error_data;
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &master_m) override {
        auto& m = static_cast<MeanComplexSign&>(master_m);
        mpi_all_reduction(result_t::value_type(N()) * result(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = (result() - m.result());
        error_data = result_t::value_type(N()) * element_wise::product(error_data, error_data);
        mpi_reduction(error_data, m.error_data);
    }
#endif
};