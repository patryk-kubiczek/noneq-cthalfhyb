#pragma once

#include "../model/cthalfhyb/Model.h"
#include "../configuration/Configuration.h"
#include "../local_weight/LocalWeightBase.h"

#include "DiscreteBath.h"


class DiscreteBathWeight : public LocalWeightBase{
public:
    DiscreteBathWeight(const Model &m, const Configuration &conf, local_weight_type type = DISCRETE_BATH) :
            LocalWeightBase(type), beta(m.beta), t_max(m.t_max), p(&m.model_params), conf(&conf),
            discrete_bath(p->eps_0_up, p->eps_up, p->U_0, p->U,
                          p->V_0_up, p->V_up, p->eps_bath_0_up, p->eps_bath_up,
                          beta, t_max, 1),
            eps_0_down(p->eps_0_down), eps_down(p->eps_down){
        generate_from_conf();
    }

    void set_n_blocks(int n_blocks);

    void generate_from_conf() override;
    cx_long_double calculate_new_trace() override;
    cx_long_double get_trace_ratio() const override;
    void finalize_update() override;
    cx_long_double get_w_loc() const override;

    arma::cx_mat density_matrix_at_t_max(int *n_down) override;

    std::ostream &print(std::ostream &os) const override;

private:

    cx_long_double trace_value = 1;
    cx_long_double new_trace_value = 1;
    int permutation_sign = 1;

    double beta;
    double t_max;
    const ModelParams *p;
    const Configuration *conf;

    DiscreteBath discrete_bath;

    double eps_0_down;
    double eps_down;
    cx_double calculate_phi(const std::vector<contour_time> &times, int n_ini) const;

    void generate_permutation_sign();

};



