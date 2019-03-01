#pragma once

#include "../Measurements.h"

class TimeStatistics : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        for(int flavor = 0; flavor < conf()->get_n_flavor(); ++flavor){
            for(const auto &t_a : conf()->get_a_times(flavor)){
                ++raw_data(t_a.branch);
            }
            for(const auto &t_c : conf()->get_c_times(flavor)){
                ++raw_data(t_c.branch);
            }
        }
    }
    void reset() override {
        if(is_master()) mean_data.zeros(3);
        else raw_data.zeros(3);
        error_data.zeros(3);
    }
    void save(const std::string &prefix, bool truncate) const override {
        const auto suffix = "times_statistics";
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + suffix + '_' + std::to_string(rank) + txt_extension,
                          result(), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        const auto suffix = "times_statistics";
        result().save(prefix + "_" + suffix + output_extension, output_file_type);
        error().save(prefix + "_" + suffix + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        using std::endl;
        os << name() << ": " << endl;
        auto statistics = result();
        auto statistics_error = error();
        for(int i = 0; i < 3; ++i){
            switch(i){
                case 0 :
                    os << " PLUS: ";
                    break;
                case 1 :
                    os << " MINUS: ";
                    break;
                case 2 :
                    os << " IMAG: ";
                    break;
            }
            os << statistics(i) << " +- " << statistics_error(i) << endl;
        }
    }
    std::string name() const override { return "Times statistics"; }
    using data_t = arma::Col<long long>;
    using result_t = arma::vec;
    result_t result() const {
        double norm = arma::sum(raw_data);
        return (is_master() ? mean_data / N() : arma::conv_to<result_t>::from(raw_data) / norm) ;
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
        auto& m = static_cast<TimeStatistics&>(global_m);
        mpi_all_reduction((N() * result()).eval(), m.mean_data);
        // Reduce squared deviation from mean
        error_data = result() - m.result();
        error_data = N() * error_data % error_data;
        mpi_reduction(error_data, m.error_data);
    }
#endif
};
