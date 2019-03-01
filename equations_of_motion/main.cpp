#include <iostream>
#include "EquationsOfMotion.h"

int main(int argc, char **argv)
{
    int grid_size = 200;
    double t_max = 5;
    int N = 10;
    // Make Hamiltonian grid twice denser than resulting grid
    double ham_dt = 0.5 * t_max / grid_size;
    // Make sure Hamiltonian t_max is at least as large as t_max in the result
    double ham_t_max = t_max + 0.1;
    // Hamiltonian should have sufficiently many elements to cover ham_t_max
    int n_ham = int(ham_t_max / ham_dt) + 1;
    // Generate random cube
    cube_t h = arma::randu<cube_t>(N, N, n_ham);
    // Ensure hermiticity
    for(int i = 0; i < n_ham; ++i){
        h.slice(i) = 0.5 * (h.slice(i) + h.slice(i).t());
    }


    EquationsOfMotion<> sol(h, ham_dt, t_max, grid_size);
    sol.run(true);
    double x = sol.check_unitarity(true);

    std::cout << x << std::endl;
}
