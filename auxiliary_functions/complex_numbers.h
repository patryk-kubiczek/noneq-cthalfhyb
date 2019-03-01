#pragma once

#include <complex>
#include "constants.h"


template <typename T>
inline double abs_error(const std::complex<T> &z, const std::complex<double> &error_z){
    using std::abs;
    return (abs(z.real()) * error_z.real() + abs(z.imag()) * error_z.imag()) / abs(z);
}

template <typename T>
inline double arg_error(const std::complex<T> &z, const std::complex<double> &error_z){
    using std::abs;
    return (abs(z.real()) * error_z.real() + abs(z.imag()) * error_z.imag()) / (abs(z) * abs(z));
}

template <typename T>
inline void print_complex_number(std::ostream &os, const std::complex<T> &z, const std::complex<double> &error_z = 0.){
    using std::abs;
    using std::arg;
    os << "(" << abs(z) << " +- " << abs_error(z, error_z) << ") * e^i("
       << arg(z) / PI << " +- " << arg_error(z, error_z) / PI << ")pi";
}
