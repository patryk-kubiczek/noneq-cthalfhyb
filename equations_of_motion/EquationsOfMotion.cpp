#include "EquationsOfMotion.h"

using namespace arma;

template <typename Hamiltonian>
void EquationsOfMotion<Hamiltonian>::run(bool verbose) {
    int grid_size = u_matrices.n_slices_1();
    for(int j = 0; j < grid_size; ++j) {
        if(verbose)
            std::cout << "Solving column " << j + 1 << "/" << grid_size << std::endl;
        if(j < grid_size - 1) {
            u_0.eye();
            solver.solve_column(j, grid_size - 1, j, u_0);
        }
        if(j > 0) {
            u_0.eye();
            solver.solve_column(j, 0, j, u_0);
        }
    }
}

template <typename Hamiltonian>
double EquationsOfMotion<Hamiltonian>::check_unitarity(bool verbose) const {
    int grid_size = u_matrices.n_slices_1();
    // Returns the largest value among the differences of determinants absolute values and unity
    double max_abs_diff = 0;
    double abs_det_diff;
    //mat errors(grid_size, grid_size);
    //cx_mat u_mat(grid_size, grid_size);

    for(int j : std::vector<int>{0, grid_size - 1}) {
        max_abs_diff = 0;
        for(int i = 0; i < grid_size; ++i) {
            abs_det_diff = std::abs(det(u_matrices.slice(i, j))) - 1;
            //errors(i, j) = abs_det_diff;
            //u_mat(i, j) = u_matrices(i,j)(0, 0);
            if(std::abs(abs_det_diff) > max_abs_diff) {
                max_abs_diff = std::abs(abs_det_diff);
            }
            if(verbose) {
                std::cout << "(i, j) = " << i << ", " << j << std::endl;
                std::cout << abs_det_diff << std::endl;
            }
        }
        std::cout << "Max error in column " << j << ": " << max_abs_diff << std::endl;
    }
    //std::cout << errors << std::endl;
    //std::cout << u_mat << std::endl;

    return max_abs_diff;
}

template <typename Hamiltonian>
void EquationsOfMotion<Hamiltonian>::save(const std::string &filename) const {
    u_matrices().save(filename);
}

