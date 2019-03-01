#pragma once

#include "contour_time.h"
#include "../auxiliary_functions/constants.h"
#include "../interpolation/Interpolation.h"

class ContourPropagator {
public:
    ContourPropagator(const cx_tetracube &U_t_t, const arma::cube &U_tau_tau,
                      const double &dt, const double &dtau, double beta)
            : U_t_t(U_t_t, dt, dt), U_tau_tau(U_tau_tau, dtau), beta(beta) {}
    arma::cx_mat operator()(const contour_time &ct_1, const contour_time &ct_2) const{
        if(ct_1 > ct_2){
            if(ct_1.branch == IMAG){
                if(ct_2.branch == IMAG)
                    return arma::conv_to<arma::cx_mat>::from(U_tau_tau(ct_1.t - ct_2.t));
                else
                    return U_tau_tau(ct_1.t) * U_t_t(0, ct_2.t);
            }
            else //s_1.branch = PLUS or MINUS
                return U_t_t(ct_1.t, ct_2.t);
        }
        else { //s_2 > s_1
            if(ct_1.branch == IMAG) //s_2.branch == IMAG
                return arma::conv_to<arma::cx_mat>::from(U_tau_tau(ct_1.t) * U_tau_tau(beta - ct_2.t));
            else { //s_1.branch = PLUS or MINUS
                if (ct_2.branch == IMAG)
                    return U_t_t(ct_1.t, 0) * U_tau_tau(beta - ct_2.t);
                else{ //s_2.branch = PLUS or MINUS
                    return U_t_t(ct_1.t, 0) * U_tau_tau(beta) * U_t_t(0, ct_2.t);
                }
            }
        }
    }
private:
    const MatrixInterpolation2D<cx_double> U_t_t;
    const MatrixInterpolation1D<double> U_tau_tau;
    const double beta;
};

class ContourPhi {
public:
    ContourPhi(const arma::cx_mat &phi_t_t, const arma::vec &phi_tau_tau,
               const double &dt, const double &dtau, double beta)
            : phi_t_t(phi_t_t, dt, dt), phi_tau_tau(phi_tau_tau, dtau), beta(beta) {}
    cx_double operator()(const contour_time &ct_1, const contour_time &ct_2) const{
        if(ct_1 > ct_2){
            if(ct_1.branch == IMAG){
                if(ct_2.branch == IMAG)
                    return phi_tau_tau(ct_1.t - ct_2.t);
                else
                    return phi_tau_tau(ct_1.t) * phi_t_t(0, ct_2.t);
            }
            else //s_1.branch = PLUS or MINUS
                return phi_t_t(ct_1.t, ct_2.t);
        }
        else { //s_2 > s_1
            //std::cout << "Bad order!" << std::endl;
            if(ct_1.branch == IMAG) //s_2.branch == IMAG
                return phi_tau_tau(ct_1.t) * phi_tau_tau(beta - ct_2.t);
            else { //s_1.branch = PLUS or MINUS
                if (ct_2.branch == IMAG)
                    return phi_t_t(ct_1.t, 0) * phi_tau_tau(beta - ct_2.t);
                else{ //s_2.branch = PLUS or MINUS
                    return phi_t_t(ct_1.t, 0) * phi_tau_tau(beta) * phi_t_t(0, ct_2.t);
                }
            }
        }
    }
private:
    const Interpolation2D<cx_double> phi_t_t;
    const Interpolation1D<double> phi_tau_tau;
    const double beta;
};
