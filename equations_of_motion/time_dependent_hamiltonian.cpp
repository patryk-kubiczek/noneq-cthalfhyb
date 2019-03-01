#include "time_dependent_hamiltonian.h"

arma::cx_mat time_dependent_hamiltonian::operator()(double t) const {
    // Simple linear interpolation
    if(t < 0)
        t = 0;
    return MatrixInterpolation1D<cx_double>{h, dt}(t);
}




