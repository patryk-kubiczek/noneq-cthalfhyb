#pragma once

#include <limits>
#include <cassert>
#include <armadillo>
#include "../auxiliary_functions/double_to_int.h"
#include "../data_types/Tetracube.h"

template <typename T>
class Interpolation1D{
public:
    Interpolation1D(const arma::Col<T> &v, const double &dt) : v(v), dt(dt) {}
    T operator()(double t) const {
        int i_lower = int_floor(t / dt);
        if(i_lower == v.n_rows - 1) i_lower = v.n_rows - 2;
        double delta_i = (t - i_lower * dt) / dt;

        assert(i_lower >= 0);
        assert(i_lower + 1 < v.n_rows);
        assert(delta_i > -1e-12 && delta_i < 1 + 1.e-12);

        return (1 - delta_i) * v(i_lower)
               + delta_i * v(i_lower + 1);
    }
private:
    const arma::Col<T> &v;
    const double &dt;
};

template <typename T>
class Interpolation2D{
public:
    Interpolation2D(const arma::Mat<T> &m, const double &dt_1, const double &dt_2) : m(m), dt_1(dt_1), dt_2(dt_2) {}
    T operator()(double t_1, double t_2) const {
        int i_lower, j_lower;
        i_lower = int_floor(t_1 / dt_1);
        if(i_lower == m.n_rows - 1) i_lower = m.n_rows - 2;
        j_lower = int_floor(t_2 / dt_2);
        if(j_lower == m.n_cols - 1) j_lower = m.n_cols - 2;
        double delta_i = (t_1 - i_lower * dt_1) / dt_1;
        double delta_j = (t_2 - j_lower * dt_2) / dt_2;

        assert(i_lower >= 0 && j_lower >= 0);
        assert(i_lower + 1 < m.n_rows && j_lower + 1 < m.n_cols);
        assert(delta_i > -1e-12 && delta_i < 1 + 1.e-12 && delta_j > -1e-12 && delta_j < 1 + 1.e-12);

        return (1 - delta_i) * (1 - delta_j) * m(i_lower, j_lower)
               + delta_i * (1 - delta_j) * m(i_lower + 1, j_lower)
               + (1 - delta_i) * delta_j * m(i_lower, j_lower + 1)
               + delta_i * delta_j * m(i_lower + 1, j_lower + 1);
    }
private:
    const arma::Mat<T> &m;
    const double &dt_1;
    const double &dt_2;
};

template <typename T>
class MatrixInterpolation1D{
public:
    MatrixInterpolation1D(const arma::Cube<T> &c, const double &dt) : c(c), dt(dt) {}
    arma::Mat<T> operator()(double t) const {
        int i_lower = int_floor(t / dt);
        if(i_lower == c.n_slices - 1) i_lower = c.n_slices - 2;
        double delta_i = (t - i_lower * dt) / dt;

        assert(i_lower >= 0);
        assert(i_lower + 1 < c.n_slices);
        assert(delta_i > -1.e-12 && delta_i < 1 + 1.e-12);

        return (1 - delta_i) * c.slice(i_lower)
               + delta_i * c.slice(i_lower + 1);
    }
private:
    const arma::Cube<T> &c;
    const double &dt;
};

template <typename T>
class MatrixInterpolation2D{
public:
    MatrixInterpolation2D(const Tetracube<T> &tc, const double &dt_1, const double &dt_2) : tc(tc), dt_1(dt_1), dt_2(dt_2) {}
    arma::Mat<T> operator()(double t_1, double t_2) const {
        int i_lower, j_lower;
        i_lower = int_floor(t_1 / dt_1);
        if(i_lower == tc.n_slices_1() - 1) i_lower = tc.n_slices_1() - 2;
        j_lower = int_floor(t_2 / dt_2);
        if(j_lower == tc.n_slices_2() - 1) j_lower = tc.n_slices_2() - 2;
        double delta_i = (t_1 - i_lower * dt_1) / dt_1;
        double delta_j = (t_2 - j_lower * dt_2) / dt_2;

        assert(i_lower >= 0 && j_lower >= 0);
        assert(i_lower + 1 < tc.n_slices_1() && j_lower + 1 < tc.n_slices_2());
        assert(delta_i > -1e-12 && delta_i < 1 + 1.e-12 && delta_j > -1e-12 && delta_j < 1 + 1.e-12);

        return (1 - delta_i) * (1 - delta_j) * tc.slice(i_lower, j_lower)
               + delta_i * (1 - delta_j) * tc.slice(i_lower + 1, j_lower)
               + (1 - delta_i) * delta_j * tc.slice(i_lower, j_lower + 1)
               + delta_i * delta_j * tc.slice(i_lower + 1, j_lower + 1);
    }
private:
    const Tetracube<T> &tc;
    const double &dt_1;
    const double &dt_2;
};