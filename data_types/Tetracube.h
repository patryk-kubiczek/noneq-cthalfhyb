#pragma once

#include <complex>
#include <string>
#include <cassert>
#include <armadillo>

template <typename T>
class Tetracube{
public:
    Tetracube() = default;
    Tetracube(const Tetracube&) = default;
    Tetracube(const arma::Cube<T> cube, int n_slices_1, int n_slices_2)
            : _n_slices_1(n_slices_1), _n_slices_2(n_slices_2), cube(cube) {}
    Tetracube(int n_rows, int n_cols, int n_slices_1, int n_slices_2)
            : _n_slices_1(n_slices_1), _n_slices_2(n_slices_2), cube(n_rows, n_cols, n_slices_1 * n_slices_2, arma::fill::zeros) {}

    void load(const std::string &name, int n_slices_1, int n_slices_2){
        _n_slices_1 = n_slices_1;
        _n_slices_2 = n_slices_2;
        cube.load(name);
    }
    void zeros(int n_rows, int n_cols, int n_slices_1, int n_slices_2) {
        _n_slices_1 = n_slices_1;
        _n_slices_2 = n_slices_2;
        cube.zeros(n_rows, n_cols, n_slices_1 * n_slices_2);
    }
    const auto& operator()() const { return cube; }
    auto& operator()(){ return cube; }
    const auto& operator()(int i, int j, int k, int l) const { return cube(i, j, k + _n_slices_1 * l); }
    auto& operator()(int i, int j, int k, int l) { return cube(i, j, k + _n_slices_1 * l); }
    const auto& slice(int k, int l) const { return cube.slice(k + l * _n_slices_1); }
    auto& slice(int k, int l) { return cube.slice(k + l * _n_slices_1); }
    int n_slices_1() const { return _n_slices_1; }
    int n_slices_2() const { return _n_slices_2; }
    Tetracube operator+(const Tetracube &other) const { return Tetracube(cube + other.cube, _n_slices_1, _n_slices_2); }
    template<typename TT>
    Tetracube operator+(const TT &number) const { return Tetracube(cube + number, _n_slices_1, _n_slices_2); }
    Tetracube operator-(const Tetracube &other) const { return Tetracube(cube - other.cube, _n_slices_1, _n_slices_2); }
    template<typename TT> Tetracube operator-(const TT &number) const{ return Tetracube(cube - number, _n_slices_1, _n_slices_2); }
    Tetracube operator%(const Tetracube &other) const { return Tetracube(cube % other.cube, _n_slices_1, _n_slices_2); }
    template<typename TT>
    Tetracube operator*(const TT &number) const{ return Tetracube(cube * number, _n_slices_1, _n_slices_2); }
    template<typename TT>
    Tetracube operator/(const TT &number) const{ return Tetracube(cube / number, _n_slices_1, _n_slices_2); }
    template<typename TT>
    void operator*=(const TT &number){ cube *= number; }
    template<typename TT>
    void operator/=(const TT &number){ cube /= number; }
private:
    int _n_slices_1 = 0;
    int _n_slices_2 = 0;
    arma::Cube<T> cube;
};

using cx_tetracube = Tetracube<std::complex<double>>;
using tetracube = Tetracube<double>;