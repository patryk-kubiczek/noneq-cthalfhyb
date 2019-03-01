#include "DiscreteBath.h"
#include "../model/contour_time.h"
#include "../gnuplot-iostream.h"
#include "../data_types/Tetracube.h"


int main(int argc, char **argv) {

    double t_max = 1;
    double beta = 1;

    double eps_0 = -3.;
    double eps = -3;
    double U_0 = 6;
    double U = 6;

    std::vector<cx_double> zeros(1, 0.);
    std::vector<cx_double> V_0 = {1.12837, 1.12837};
    std::vector<cx_double> V = {1.12837, 1.12837};
    std::vector<double> eps_bath_0 = {0, 0};
    std::vector<double> eps_bath = {-1, 1};

    int n_ini = 0;
    std::vector<contour_time> times = {{PLUS, 0.1 * t_max},{PLUS, 0.5 * t_max},{IMAG, 0.5 * beta},{IMAG, 0.9 * beta}};

    int n_block = 1;

    DiscreteBath discrete_bath(eps_0, eps, U_0, U, V_0, V, eps_bath_0, eps_bath, beta, t_max, n_block);



    discrete_bath.generate_u(times, n_ini);
    discrete_bath.generate_QDT();

    std::cout << std::exp(discrete_bath.log_z()) << ", " << std::exp(discrete_bath.log_z_without_QDT()) << std::endl;
    std::cout << discrete_bath.density_matrix_at_t_max() << std::endl << discrete_bath.density_matrix_at_t_max_without_QDT() << std::endl;

}