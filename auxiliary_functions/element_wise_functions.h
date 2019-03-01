#pragma once

#include <complex>
#include "../data_types/Tetracube.h"

namespace element_wise {
    template <typename T>
    inline std::complex<T> sqrt(const std::complex<T> &z){
        return {std::sqrt(z.real()), std::sqrt(z.imag())};
    }

    template <typename T>
    inline std::complex<T> product(const std::complex<T> &lhs, const std::complex<T> &rhs){
        return {lhs.real() * rhs.real(), lhs.imag() * rhs.imag()};
    }

    template <typename T>
    inline arma::Col<std::complex<T>> product(const arma::Col<std::complex<T>> &lhs,
                                              const arma::Col<std::complex<T>> &rhs){
        return arma::Col<std::complex<T>>{arma::real(lhs) % arma::real(rhs),
                                          arma::imag(lhs) % arma::imag(rhs)};
    }

    template <typename T>
    inline arma::Mat<std::complex<T>> product(const arma::Mat<std::complex<T>> &lhs,
                                              const arma::Mat<std::complex<T>> &rhs){
        return arma::Mat<std::complex<T>>{arma::real(lhs) % arma::real(rhs),
                                          arma::imag(lhs) % arma::imag(rhs)};
    }

    template <typename T>
    inline arma::Cube<std::complex<T>> sqrt(const arma::Cube<std::complex<T>> &c){
        using arma::sqrt;
        return arma::Cube<std::complex<T>>{(sqrt(arma::real(c))), sqrt(arma::imag(c))};
    }

    template <typename T>
    inline Tetracube<std::complex<T>> sqrt(const Tetracube<std::complex<T>> &tc){
        using arma::sqrt;
        return {arma::Cube<std::complex<T>>{sqrt(arma::real(tc())), sqrt(arma::imag(tc()))},
                tc.n_slices_1(), tc.n_slices_2()};
    }

    template <typename T>
    inline arma::Cube<std::complex<T>> product(const arma::Cube<std::complex<T>> &lhs,
                                                      const arma::Cube<std::complex<T>> &rhs){
        return arma::Cube<std::complex<T>>{arma::real(lhs) % arma::real(rhs),
                                            arma::imag(lhs) % arma::imag(rhs)};
    }

    template <typename T>
    inline Tetracube<std::complex<T>> product(const Tetracube<std::complex<T>> &lhs,
                                                      const Tetracube<std::complex<T>> &rhs){
        return {arma::Cube<std::complex<T>>{arma::real(lhs()) % arma::real(rhs()),
                                            arma::imag(lhs()) % arma::imag(rhs())},
                lhs.n_slices_1(), lhs.n_slices_2()};
    }
}


