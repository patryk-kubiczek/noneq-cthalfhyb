#pragma once

#include "../Measurements.h"

// This measurement is valid only for CT-1/2-HYB

class SpinUp : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        int n_down;
        arma::cx_mat contribution = -const_cast<local_weight_t*>(local_weight())->density_matrix_at_t_max(&n_down);
        contribution.diag() += 1;
        // Density matrix elements
        raw_data[0] += contribution * conf()->get_sign();
        // Double occupancy
        raw_data[1](0) += double(n_down) * contribution(0) * conf()->get_sign();
        DEBUG( std::cout << "Measuring spin up quantities"  << std::endl;
        std::cout << "n_down = " << n_down  << std::endl;
        std::cout << "n_up = " << contribution[0] << std::endl; )
    }
    void reset() override {
        if(is_master()) {
            mean_data.clear();
            mean_data.emplace_back(arma::cx_mat(n_bath() + 1, n_bath() + 1, arma::fill::zeros));
            mean_data.emplace_back(arma::cx_mat(1, 1, arma::fill::zeros));
        }
        else {
            raw_data.clear();
            raw_data.emplace_back(arma::cx_mat(n_bath() + 1, n_bath() + 1, arma::fill::zeros));
            raw_data.emplace_back(arma::cx_mat(1, 1, arma::fill::zeros));
        }
        error_data.clear();
        error_data.emplace_back(arma::cx_mat(n_bath() + 1, n_bath() + 1, arma::fill::zeros));
        error_data.emplace_back(arma::cx_mat(1, 1, arma::fill::zeros));
    }
    void save(const std::string &prefix, bool truncate) const override {
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        save_to_row::save(prefix + "_" + "spin_up_density_matrix_elements" + '_' + std::to_string(rank) + txt_extension,
                          result(0), {mean_sign().real(), double(N()), double(rank)}, truncate);
        save_to_row::save(prefix + "_" + "double_occupancy" + '_' + std::to_string(rank) + txt_extension,
                          result(1)(0), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        result(0).save(prefix + "_" + "spin_up_density_matrix_elements" + output_extension, output_file_type);
        error(0).save(prefix + "_" + "spin_up_density_matrix_elements" + "_error" + output_extension, output_file_type);
        result(1).save(prefix + "_" + "double_occupancy" + output_extension, output_file_type);
        error(1).save(prefix + "_" + "double_occupancy" + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        os << name() << ":" << std::endl;
        os << " n_up(t_max) = " << std::real(result(0)(0)) << " +- " << std::real(error(0)(0)) << std::endl;
        os << " d(t_max) = " << std::real(result(1)(0)) << " +- " << std::real(error(1)(0)) << std::endl;
    }
    std::string name() const override { return "Spin-up impurity-bath density matrix elements at t_max"; }
    using data_t = std::vector<arma::cx_mat>;
    using result_t = data_t;
    result_t::value_type result(int index) const {
        return (is_master() ? mean_data[index] : raw_data[index] / mean_sign()) /  N();
    }
    result_t::value_type error(int index) const {
        return arma::sqrt(error_data[index] / N());
    }
private:
    data_t raw_data;
    result_t mean_data;
    result_t error_data;
    double dt;
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &global_m) override {
        auto& m = static_cast<SpinUp&>(global_m);
        for(int i = 0; i < 2; ++i){
            mpi_all_reduction((N() * result(i)).eval(), m.mean_data[i]);
            // Reduce squared deviation from mean
            error_data[i] = result(i) - m.result(i);
            error_data[i] = element_wise::product(error_data[i], error_data[i]) * N();
            mpi_reduction(error_data[i], m.error_data[i]);
        }
    }
#endif
};
