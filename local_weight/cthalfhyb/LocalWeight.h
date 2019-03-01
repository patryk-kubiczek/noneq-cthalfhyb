#pragma once

#include <list>
#include <functional>
#include <memory>
#include <algorithm>
#include <ostream>
#include <bitset>

#include "../../model/ContourPropagatorWrapper.h"
#include "../../configuration/Configuration.h"
#include "../../auxiliary_functions/get_permutation_sign.h"
#include "../../auxiliary_functions/container_functions.h"
#include "../LocalWeightBase.h"

class LocalWeight : public LocalWeightBase {
public:
    LocalWeight(const Model &m, const Configuration &conf, local_weight_type type);

    void generate_from_conf() override;
    cx_long_double calculate_new_trace() override;
    cx_long_double get_trace_ratio() const override;
    // Run this after all the confirmations if move is accepted
    void finalize_update() override;
    // Run this if move is not accepted
    void cancel_operations();
    cx_long_double get_w_loc() const override;
    arma::cx_mat density_matrix_at_t_max(int *n_down_at_t_max) override;

    std::ostream &print(std::ostream &os) const override;

    class EvolutionMatrix {
    public:
        EvolutionMatrix(const contour_time& end_time, const contour_time& start_time, arma::cx_mat u_matrix, cx_double phi)
                : start_time(start_time), end_time(end_time), u_matrix(std::move(u_matrix)), phi(phi) {}
        const arma::cx_mat& get() const { return u_matrix; }
        cx_double get_phi() const { return phi; }
        const contour_time& get_start_time() const { return start_time; }
        const contour_time& get_end_time() const { return end_time; }
    private:
        contour_time end_time;
        contour_time start_time;
        arma::cx_mat u_matrix;
        cx_double phi;
    };

    using evolution_matrices_t = std::list<EvolutionMatrix>;
    using matrix_iters_t = std::vector<evolution_matrices_t::iterator>;

    // relative index is the index of d^dagger in case n_down = 0 and of d in case n_down = 1
    // n_down always corresponds to the evolution between the pair operators, either removed or added

    void attempt_add_pair(int n_down, const contour_time &end_time, const contour_time &start_time);
    void attempt_remove_pair(int n_down, int relative_index);

    void attempt_add_wound_pair(int n_down, const contour_time &end_time, const contour_time &start_time);
    void attempt_remove_wound_pair(int n_down);

    void attempt_switch_n(int new_n_down);

    void confirm_add_pair(int n_down, int c_relative_index, int a_relative_index);
    void confirm_remove_pair(int n_down, int c_relative_index, int a_relative_index);

    void confirm_add_wound_pair(int n_down);
    void confirm_remove_wound_pair(int n_down);

private:
    double t_max;
    double beta;

    ContourUWrapper<0> u_0;
    ContourUWrapper<1> u_1;
    ContourPhi1Wrapper phi_1;

    const Configuration *conf;

    cx_long_double trace_value = 1;
    cx_long_double new_trace_value = 1;
    int permutation_sign = 1;

    EvolutionMatrix make_evolution_matrix(int n_down, const contour_time &end_time, const contour_time &start_time) const;

    evolution_matrices_t evolution_matrices;

    // Those vectors store iterators to evolution matrices PRECEDING a given operator
    matrix_iters_t c_matrix_iters;
    matrix_iters_t a_matrix_iters;

    std::vector<evolution_matrices_t> old_evolution_matrices;
    std::vector<std::pair<evolution_matrices_t::iterator, evolution_matrices_t::iterator>> new_mat_iters;

    evolution_matrices_t::iterator temporary_iter_removal;
    evolution_matrices_t::iterator temporary_iter_addition;

    int get_current_order() const;
    void generate_permutation_sign();
    cx_long_double determinant() const;
};




