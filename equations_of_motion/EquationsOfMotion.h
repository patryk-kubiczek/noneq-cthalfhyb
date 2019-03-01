#pragma once

#include "Solver.h"

template <typename Hamiltonian = time_dependent_hamiltonian>
class EquationsOfMotion {
    using ham_t = Hamiltonian;
public:
    EquationsOfMotion(const std::string &ham_filename, double ham_dt, double t_max, int grid_size)
        : ham(ham_filename, ham_dt),
          u_matrices(ham(0).n_rows, ham(0).n_cols, grid_size, grid_size),
          solver(grid_size, t_max, ham, u_matrices),
          u_0(arma::size(ham(0))) {
        assert(ham.t_max() >= t_max);
    }
    EquationsOfMotion(const arma::cx_cube &h, double ham_dt, double t_max, int grid_size)
        : ham(h, ham_dt),
          u_matrices(h.n_rows, h.n_cols, grid_size, grid_size),
          solver(grid_size, t_max, ham, u_matrices),
          u_0(arma::size(h.slice(0))) {
        assert(ham.t_max() >= t_max);
    }
    void run(bool verbose = false);
    double check_unitarity(bool verbose = false) const;
    void save(const std::string &result_filename) const;
    const cx_tetracube &get_result() const {
        return u_matrices;
    }
private:
    const ham_t ham;
    cx_tetracube u_matrices;
    Solver<ham_t> solver;
    arma::cx_mat u_0;
};

template class EquationsOfMotion<time_dependent_hamiltonian>;






