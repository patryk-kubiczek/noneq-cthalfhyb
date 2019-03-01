#pragma once

#include <ostream>
#include "../auxiliary_functions/constants.h"
#include "../configuration/Configuration.h"
#include "../auxiliary_functions/matrix_functions.h"
#include "../model/analytic_expressions.h"
#include "../model/ContourDeltaWrapper.h"

class BathWeight {
public:
    BathWeight(const Model &model, const Configuration &conf, bath_weight_type type)
            : delta(ContourDeltaWrapper::create(model, type)), conf(&conf) {}

    cx_double get_w_bath();
    cx_double get_det_ratio();
    void generate_from_conf();

    void attempt_insertion(int flavor, const contour_time &t_c, const contour_time &t_a);
    void attempt_removal(int i, int j);
    void attempt_shift_t_a(int flavor, int j, const contour_time &new_t_a);
    void attempt_shift_t_c(int flavor, int i, const contour_time &new_t_c);

    void confirm_insertion(int i, int j);
    void confirm_removal(int i, int j);
    void confirm_shift_t_a(int flavor, int j, int winding_direction);
    void confirm_shift_t_c(int flavor, int j, int winding_direction);

    cx_double cut_hybridization_line(int i, int j) const;

    friend std::ostream &operator<<(std::ostream &os, const BathWeight &weight);

private:
    ContourDeltaWrapper delta;
    const Configuration *conf;

    cx_double det_value;
    arma::cx_mat M;

    arma::cx_colvec col;
    arma::cx_rowvec row;
    arma::Mat<cx_double>::fixed<1, 1> scalar;

    cx_double det_ratio;

    cx_double get_sign_factor() const;
};


