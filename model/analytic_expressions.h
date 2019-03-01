#pragma once

#include <armadillo>
#include <cmath>
#include "../auxiliary_functions/constants.h"
#include "contour_time.h"
#include "params.h"

class AnalyticDelta{
public:
    AnalyticDelta(double beta, const ModelParams &m, int spin) : beta(beta), m(&m), spin(spin) {
        //assert(m->eps_bath_down.size() == 0);
    }
    cx_double operator()(const contour_time &ct_1, const contour_time &ct_2) const {
        double eps_bath_0 = (spin == 0 ? m->eps_bath_0_up[0] : m->eps_bath_0_down[0]);
        double eps_bath = (spin == 0 ? m->eps_bath_up[0] : m->eps_bath_down[0]);
        cx_double V_0 = (spin == 0 ? m->V_0_up[0] : m->V_0_down[0]);
        cx_double V = (spin == 0 ? m->V_up[0] : m->V_down[0]);

        double fermi = 1. / (1 + std::exp(beta * eps_bath_0));
        cx_double time_diff = ct_1.t * (ct_1.branch == IMAG ? -I : 1.)
                              - ct_2.t * (ct_2.branch == IMAG ? -I : 1.);
        cx_double exponential = std::exp(-I * time_diff.real() * eps_bath + time_diff.imag() * eps_bath_0);
        return (ct_1.branch == IMAG ? V_0 : V)
               * -I *((ct_1 > ct_2 ? 1 : 0) - fermi) * exponential
               * (ct_2.branch == IMAG ? V_0 : V);
    }

private:
    int spin;
    double beta;
    const ModelParams *m;

};

template <unsigned n>
class AnalyticU{
public:
    explicit AnalyticU(double beta, const ModelParams &m) : beta(beta), m(&m) {
        //assert(m->eps_bath_down.size() == 0);
        //assert(m->eps_bath_down[0] == 0 && m->eps_bath_0_down[0] == 0);
    }
    arma::cx_mat operator()(const contour_time &ct_1, const contour_time &ct_2) const {
        if(ct_1 < ct_2)
            return (*this)(ct_1, {PLUS, 0}) * (*this)({IMAG, beta}, ct_2);
//        arma::mat h_0 = { {m->eps_0_up + n * m->U_0, std::abs(m->V_0_down[0])},
//                          {std::abs(m->V_0_down[0]), m->eps_bath_0_down[0]}      };
//        arma::mat h = { {m->eps_up + n * m->U, std::abs(m->V_down[0])},
//                          {std::abs(m->V_down[0]), m->eps_bath_down[0]}      };
        cx_double time_diff = ct_1.t * (ct_1.branch == IMAG ? -I : 1.)
                              - ct_2.t * (ct_2.branch == IMAG ? -I : 1.);
        double t = time_diff.real();
        double tau = time_diff.imag();
//        arma::cx_mat res_0 = arma::conv_to<arma::cx_mat>::from(arma::expmat_sym(tau * h_0))
//               * arma::expmat(-I * t * h);
        //return res_0;
        using std::exp; using std::sqrt;
        double deps_0 = m->eps_0_down + n * m->U_0;
        double deps = m->eps_down + n * m->U;
        double eeps_0 = std::sqrt(deps_0 * deps_0 + 4. * std::abs(m->V_0_down[0]) * std::abs(m->V_0_down[0]));
        double eeps = std::sqrt(deps * deps + 4. * std::abs(m->V_down[0]) * std::abs(m->V_down[0]));
        arma::cx_mat res(2, 2, arma::fill::eye);
        if(tau != 0){
            using std::sinh; using std::cosh;
            res *= exp(deps_0 * tau / 2.);
            double a = cosh(eeps_0 * tau / 2.);
            double b = sinh(eeps_0 * tau / 2.);
            res *= {{a + deps_0 / eeps_0 * b, 2. * m->V_0_down[0] / eeps_0 * b},
                    {2. * m->V_0_down[0] / eeps_0 * b, a - deps_0 / eeps_0 * b}};
        }
        if(t != 0){
            using std::cos; using std::sin;
            res *= exp(-I * deps * t / 2.);
            double a = cos(eeps * t / 2.);
            double b = sin(eeps * t / 2.);
            res *= {{a - deps / eeps * b * I, -2. * m->V_down[0] / eeps * b * I},
                    {-2. * m->V_down[0] / eeps * b * I, a + deps / eeps * b * I}};
        }
        return res;
    }
private:
    double beta;
    const ModelParams *m;
};


class AnalyticPhi1{
public:
    explicit AnalyticPhi1(double beta, const ModelParams &m) : beta(beta), m(&m) {}
    cx_double operator()(const contour_time &ct_1, const contour_time &ct_2) const {
        if(ct_1 < ct_2)
            return (*this)(ct_1, {PLUS, 0}) * (*this)({IMAG, beta}, ct_2);
        cx_double time_diff = ct_1.t * (ct_1.branch == IMAG ? -I : 1.)
                              - ct_2.t * (ct_2.branch == IMAG ? -I : 1.);
        return std::exp(-I * time_diff.real() * m->eps_up + time_diff.imag() * m->eps_0_up);
    }
private:
    double beta;
    const ModelParams *m;
};

class AnalyticP{
public:
    explicit AnalyticP(double beta, const ModelParams &m) : beta(beta), m(&m) {}
    arma::cx_mat operator()(const contour_time &ct_1, const contour_time &ct_2) const {
        if(ct_1 < ct_2)
            return (*this)(ct_1, {PLUS, 0}) * (*this)({IMAG, beta}, ct_2);
        arma::vec energies_0 = {0, m->eps_0_up, m->eps_0_down, m->U_0 + m->eps_0_up + m->eps_0_down};
        arma::vec energies = {0, m->eps_up, m->eps_down, m->U + m->eps_up + m->eps_down};
        cx_double time_diff = ct_1.t * (ct_1.branch == IMAG ? -I : 1.)
                              - ct_2.t * (ct_2.branch == IMAG ? -I : 1.);
        arma::cx_mat result(4, 4, arma::fill::zeros);
        for(int i = 0; i < 4; ++i){
            result(i, i) = std::exp(-I * time_diff.real() * energies(i) + time_diff.imag() * energies_0(i));
        }
        return result;
    }

private:
    double beta;
    const ModelParams *m;
};

