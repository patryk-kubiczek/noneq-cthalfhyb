#pragma once

#include <armadillo>
#include <cmath>
#include <string>
#include <cassert>

#include "../auxiliary_functions/constants.h"
#include "../interpolation/Interpolation.h"

using cube_t = arma::cx_cube;

class time_dependent_hamiltonian{
public:
    time_dependent_hamiltonian(const arma::cx_cube &h, double dt) : h(h), dt(dt) {}
    time_dependent_hamiltonian(const std::string &filename, double dt) : h(load_h(filename)), dt(dt) {}
    inline double t_max() const{
        return h.n_slices * dt;
    };
    arma::cx_mat operator()(double t) const;
private:
    const arma::cx_cube h;
    const double dt;
    static arma::cx_cube load_h(const std::string &filename){
        arma::cx_cube h_tmp;
        h_tmp.load(filename, arma::arma_binary);
        return h_tmp;
    }
};



