#pragma once

#include <cmath>
#include <complex>

#ifndef NDEBUG
    #define DEBUG(x) x
#else
    #define DEBUG(x)
#endif

using cx_double = std::complex<double>;
using cx_long_double = std::complex<long double>;

constexpr cx_double I{0, 1};
constexpr double PI{std::atan(1) * 4};

