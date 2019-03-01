#pragma once

#include <armadillo>
#include "../auxiliary_functions/constants.h"

enum local_weight_type { NUMERIC_LW, ANALYTIC_LW, DISCRETE_BATH };

class LocalWeightBase {
public:
    explicit LocalWeightBase(local_weight_type type) : type(type) {}
    virtual void generate_from_conf() = 0;
    virtual cx_long_double calculate_new_trace() = 0;
    virtual cx_long_double get_trace_ratio() const = 0;
    virtual void finalize_update() = 0;
    virtual cx_long_double get_w_loc() const = 0;
#ifdef CTHALFHYB_QMC
    virtual arma::cx_mat density_matrix_at_t_max(int *n_down_at_t_max) = 0;
#endif
    virtual std::ostream &print(std::ostream &os) const = 0;
    friend std::ostream &operator<<(std::ostream &os, const LocalWeightBase &weight) { return weight.print(os); };
    const local_weight_type type;
};

