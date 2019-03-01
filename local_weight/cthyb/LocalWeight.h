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

    std::ostream &print(std::ostream &os) const override;

    // Special members for Local Weight
    class LocalMatrix {
    public:
        virtual const arma::cx_mat& get() const = 0;
        virtual const contour_time& get_start_time() const = 0;
        virtual const contour_time& get_end_time() const = 0;
        virtual ~LocalMatrix() = default;
    };

    using local_matrices_t = std::list<std::unique_ptr<LocalMatrix>>;
    using matrix_iters_t = std::vector<local_matrices_t::iterator>;

    void attempt_add_operator(operator_type type, int flavor, const contour_time &time);
    void attempt_remove_operator(operator_type op_type, int flavor, int relative_index);

    void confirm_add_operator(operator_type op_type, int flavor, int relative_index);
    void confirm_remove_operator(operator_type op_type, int flavor, int relative_index);

private:
    double t_max;
    double beta;

    ContourPropagatorWrapper p;  // bare_propagator
    const arma::cx_cube *c_operators;
    const arma::cx_cube *a_operators;

    const Configuration *conf;

    cx_long_double trace_value = 1;
    cx_long_double new_trace_value = 1;
    int permutation_sign = 1;

    class LocalFermionicOperator : public LocalMatrix {
    public:
        LocalFermionicOperator(const contour_time& time, const arma::cx_mat& op_matrix)
                : time(time), op_matrix(op_matrix) {}
        const arma::cx_mat& get() const override { return op_matrix; }
        const contour_time& get_start_time() const override { return time; }
        const contour_time& get_end_time() const override { return time; }

    private:
        contour_time time;
        const arma::cx_mat& op_matrix;

    };

    class LocalPropagator : public LocalMatrix {
    public:
        LocalPropagator(const contour_time& end_time, const contour_time& start_time, arma::cx_mat p_matrix)
                : start_time(start_time), end_time(end_time), p_matrix(std::move(p_matrix)) {}
        const arma::cx_mat& get() const override { return p_matrix; }
        const contour_time& get_start_time() const override { return start_time; }
        const contour_time& get_end_time() const override { return end_time; }
    private:
        contour_time start_time;
        contour_time end_time;
        arma::cx_mat p_matrix;
    };

    std::unique_ptr<LocalFermionicOperator> make_local_operator(operator_type type, int flavor, const contour_time& time);
    std::unique_ptr<LocalPropagator> make_local_propagator(const contour_time &end_time, const contour_time &start_time);

    local_matrices_t local_matrices;

    std::vector<matrix_iters_t> c_matrix_iters;
    std::vector<matrix_iters_t> a_matrix_iters;

    std::vector<local_matrices_t> old_local_matrices;
    std::vector<std::pair<local_matrices_t::iterator, local_matrices_t::iterator>> new_mat_iters;

    local_matrices_t::iterator temporary_c_iter;
    local_matrices_t::iterator temporary_a_iter;

    bool is_zero_order_line_present() const;
    int get_contributing_matrix_element() const;
    void generate_permutation_sign();
    cx_double trace_of_propagators(int matrix_element);
};

