#pragma once

#include "../Measurements.h"

class MatsubaraGreenFunctions : public Measurements::MeasurementBase {
public:
    using MeasurementBase::MeasurementBase;

    void measure() override {
        cx_double contribution;
        int bin_a, bin_c;
        int i = 0;
        int j = 0;

#ifdef CTHALFHYB_QMC
        bool operators_on_real_branch = conf()->has_operators_on_real_branch();
        operator_type max_op; int max_flavor; int max_position;
        if(operators_on_real_branch)
            conf()->maximum_real_time(max_op, max_flavor, max_position);
#endif

        for(int flavor_c = 0; flavor_c < conf()->get_n_flavor(); flavor_c++){
            for(int pos_c = 0; pos_c < conf()->get_flavor_order(flavor_c); ++pos_c){
                j = 0;
                for(int flavor_a = 0; flavor_a < conf()->get_n_flavor(); flavor_a++){
                    for(int pos_a = 0; pos_a < conf()->get_flavor_order(flavor_a); ++pos_a){
#ifdef CTHALFHYB_QMC
                        if(operators_on_real_branch) {
                            if(max_op == creation){
                                if(flavor_c != max_flavor || pos_c != max_position) { ++j; continue; }
                            }
                            else if (max_op == annihilation) {
                                if(flavor_a != max_flavor || pos_a != max_position) { ++j; continue; }
                            }
                        }
#endif

                        const auto &t_c = conf()->get_contour_time(creation, flavor_c, pos_c);
                        const auto &t_a = conf()->get_contour_time(annihilation, flavor_a, pos_a);

                        contribution = bath_weight()->cut_hybridization_line(i, j);
                        contribution /= -t_c.integration_sign() * t_a.integration_sign();
                        contribution *= conf()->get_sign();

                        DEBUG( std::cout << contribution << std::endl;
                        std::cout << "branch " << t_a.branch << ": " << t_a.t << std:: endl;
                        std::cout << "branch " << t_c.branch << ": " << t_c.t << std:: endl; )

                        if(t_a.branch == IMAG || t_c.branch == IMAG){
                            // tau_t
                            if(t_c.branch != IMAG){
                                DEBUG( std::cout << "Г" << std::endl; )
                                int bin_c = int_round(t_c.t * (n_times() - 1) / t_max());
                                raw_data[0].slice(flavor_a, flavor_c).col(bin_c)
                                        += contribution * arma::exp(I * wn * t_a.t);
                            }
                            else if(t_a.branch != IMAG){
                                DEBUG( std::cout << "Г conjugate" << std::endl; )
                                int bin_a = int_round(t_a.t * (n_times() - 1) / t_max());
                                raw_data[0].slice(flavor_c, flavor_a).col(bin_a)
                                        += -std::conj(contribution) * arma::exp(-I * wn * t_c.t);
                            }
                                // tau
                            else{
                                DEBUG( std::cout << "M" << std::endl; )
                                // Both are IMAG times
                                raw_data[1].slice(flavor_a, flavor_c).col(0)
                                        += contribution * arma::exp(I * wn * (t_a.t - t_c.t));
                            }
                        }
                        ++j;
                    }
                }
                ++i;
            }
        }
    }
    void reset() override {
        auto reset_data = [this](data_t &data){
            int n_flavor = conf()->get_n_flavor();
            data.clear();
            data.emplace_back(cx_tetracube{n_matsubara(), n_times(), n_flavor, n_flavor});
            data.emplace_back(cx_tetracube{n_matsubara(), 1, n_flavor, n_flavor});
        };
        if(is_master()) {
            reset_data(mean_data);
        }
        else{
            reset_data(raw_data);
        }
        reset_data(error_data);

        wn.resize(n_matsubara());
        for(int n = 0; n < n_matsubara(); ++n){
            wn(n) = (2 * n + 1) * PI / beta();
        }
    }
    void save(const std::string &prefix, bool truncate) const override {
        auto suffix = [](const std::string &gf_type){ return "gf_" + gf_type; };
        int rank = 0;
#ifdef USE_MPI
        rank = mpi_rank();
#endif
        if(t_max() > 0){
            save_to_row::save(prefix + "_" + suffix("iwn_t") + '_' + std::to_string(rank) + txt_extension,
                              result(0), {mean_sign().real(), double(N()), double(rank)}, truncate);
        }
        save_to_row::save(prefix + "_" + suffix("iwn") + '_' + std::to_string(rank) + txt_extension,
                          result(1), {mean_sign().real(), double(N()), double(rank)}, truncate);
    }
    void master_save(const std::string &prefix) const override {
        auto suffix = [](const std::string &gf_type){ return "gf_" + gf_type; };
        if(t_max() > 0){
            result(0)().save(prefix + "_" + suffix("iwn_t") + output_extension, output_file_type);
            error(0)().save(prefix + "_" + suffix("iwn_t") + "_error" + output_extension, output_file_type);
        }
        result(1)().save(prefix + "_" + suffix("iwn") + output_extension, output_file_type);
        error(1)().save(prefix + "_" + suffix("iwn") + "_error" + output_extension, output_file_type);
    }
    void print_normalized(std::ostream &os) const override {
        os << "Matsubara Green functions have been measured." << std::endl;
    }
    std::string name() const override { return "Matsubara Green functions: Г iwn, I iwn"; }
    using data_t = std::vector<cx_tetracube>; // vector of 2-dim containers of GF matrices in time space
    using result_t = data_t;
    result_t::value_type result(int index) const {
        auto res = (is_master() ? mean_data[index] : raw_data[index] / mean_sign()) /  N();
        // Correct for the smaller bins at the edges of the time domain
        if(!is_master()){
            for(int a = 0; a < conf()->get_n_flavor(); ++a){
                for(int b = 0; b < conf()->get_n_flavor(); ++b){
                    int n_cols = res.slice(a, b).n_cols;
                    if(n_cols > 1){
                        res.slice(a, b).col(0) *= 2;
                        res.slice(a, b).col(n_cols - 1) *= 2;
                    }
                }
            }
            if(index == 0){
                res *=  (n_times() - 1) / t_max() / 4;
            }
            else if(index == 1) {
                res *= 1. / beta();
            }
        }
        return res;

    }
    result_t::value_type error(int i) const {
        return element_wise::sqrt(error_data[i] / N());
    }
private:
    data_t raw_data;
    result_t mean_data;
    result_t error_data;

    arma::vec wn;
public:
#ifdef USE_MPI
    void reduce(MeasurementBase &global_m) override {
        auto& m = static_cast<MatsubaraGreenFunctions&>(global_m);
        for(int i = 0; i < 2; ++i){
            if(t_max() == 0 && i < 1) continue;
            mpi_all_reduction(result(i) * N(), m.mean_data[i]);
            // Reduce squared deviation from mean
            error_data[i] = result(i) - m.result(i);
            error_data[i] = element_wise::product(error_data[i], error_data[i]) * N();
            mpi_reduction(error_data[i], m.error_data[i]);
        }
    }
#endif
};

