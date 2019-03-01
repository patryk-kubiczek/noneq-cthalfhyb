#pragma once

#include <armadillo>
#include "../auxiliary_functions/constants.h"
#include "../interpolation/Interpolation.h"
#include "contour_time.h"

class ContourDelta {
public:
    ContourDelta(const arma::field<arma::cx_mat> &delta_greater_t_t_data,
              const arma::field<arma::cx_mat> &delta_lesser_t_t_data,
              const arma::field<arma::cx_mat> &delta_tau_t_data,
              const double &dt, const double &dtau, double beta)
            : delta_greater_t_t_data(&delta_greater_t_t_data),
              delta_lesser_t_t_data(&delta_lesser_t_t_data),
              delta_tau_t_data(&delta_tau_t_data),
              beta(beta), dt(&dt), dtau(&dtau) {}
    cx_double operator()(int a, int b, const contour_time &ct_1, const contour_time &ct_2) const {
        if(ct_1 > ct_2) {
            if(ct_1.branch == IMAG) {
                if(ct_2.branch == IMAG)
                    return delta_tau_t(a, b)(ct_1.t - ct_2.t, 0);
                else
                    return delta_tau_t(a, b)(ct_1.t, ct_2.t);
            }
            else // ct_1.branch = PLUS or MINUS
                return delta_greater_t_t(a, b)(ct_1.t, ct_2.t);
        }
        else{
            if(ct_1.branch == IMAG) //ct_2.branch == IMAG
                return -delta_tau_t(a, b)(beta  + ct_1.t - ct_2.t, 0);
            else { //ct_1.branch = PLUS or MINUS
                if (ct_2.branch == IMAG)
                    return std::conj(delta_tau_t(b, a)(beta - ct_2.t, ct_1.t));
                else //ct_2.branch = PLUS or MINUS
                    return delta_lesser_t_t(a, b)(ct_1.t, ct_2.t);
            }
        }
    }

private:
    const arma::field<arma::cx_mat> *delta_greater_t_t_data;
    const arma::field<arma::cx_mat> *delta_lesser_t_t_data;
    const arma::field<arma::cx_mat> *delta_tau_t_data;
    Interpolation2D<cx_double> delta_greater_t_t(int a, int b) const {
        return {(*delta_greater_t_t_data)(a, b), *dt, *dt};
    }
    Interpolation2D<cx_double> delta_lesser_t_t(int a, int b) const {
        return {(*delta_lesser_t_t_data)(a, b), *dt, *dt};
    }
    Interpolation2D<cx_double> delta_tau_t(int a, int b) const {
        return {(*delta_tau_t_data)(a, b), *dtau, *dt};
    }
    double beta;
    const double *dt;
    const double *dtau;
};

// MORE EFFICIENT IMPLEMENTATION FOR CT-1/2-HYB

class ContourDeltaSingle {
public:
    ContourDeltaSingle(const arma::cx_mat &delta_greater_t_t,
                 const arma::cx_mat &delta_lesser_t_t,
                 const arma::cx_mat &delta_tau_t,
                 const double &dt, const double &dtau, double beta)
            : delta_greater_t_t(delta_greater_t_t, dt, dt),
              delta_lesser_t_t(delta_lesser_t_t, dt, dt),
              delta_tau_t(delta_tau_t, dtau, dt),
              beta(beta) {}
    cx_double operator()(const contour_time &ct_1, const contour_time &ct_2) const {
        if(ct_1 > ct_2) {
            if(ct_1.branch == IMAG) {
                if(ct_2.branch == IMAG)
                    return delta_tau_t(ct_1.t - ct_2.t, 0);
                else
                    return delta_tau_t(ct_1.t, ct_2.t);
            }
            else // ct_1.branch = PLUS or MINUS
                return delta_greater_t_t(ct_1.t, ct_2.t);
        }
        else{
            if(ct_1.branch == IMAG) //ct_2.branch == IMAG
                return -delta_tau_t(beta  + ct_1.t - ct_2.t, 0);
            else { //ct_1.branch = PLUS or MINUS
                if (ct_2.branch == IMAG)
                    return std::conj(delta_tau_t(beta - ct_2.t, ct_1.t));
                else //ct_2.branch = PLUS or MINUS
                    return delta_lesser_t_t(ct_1.t, ct_2.t);
            }
        }
    };
    cx_double operator()(int a, int b, const contour_time &ct_1, const contour_time &ct_2) const {
        return operator()(ct_1, ct_2);
    };
private:
    Interpolation2D<cx_double> delta_greater_t_t;
    Interpolation2D<cx_double> delta_lesser_t_t;
    Interpolation2D<cx_double> delta_tau_t;
    double beta;
};

