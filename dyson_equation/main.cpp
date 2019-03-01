#include "DysonEquation.h"
#include "../noneq_qmc/Model/contour_time.h"
#include "../gnuplot-iostream.h"


int main(int argc, char **argv)
{

    double t_max = 1;
    double beta = 10;

    double eps_0 = 2.;
    double eps = -2;
    double U_0 = 4;
    double U = 4;

    std::vector<complex_t> zeros(1, 0.);
    std::vector<complex_t> V_0 = {1.};
    std::vector<complex_t> V = {1.};
    std::vector<double> eps_bath_0 = {-5};
    std::vector<double> eps_bath = {-4};

    int n_ini = 0;
    std::vector<contour_time> times = {{PLUS, 0.1 * t_max},{PLUS, 0.5 * t_max},{IMAG, 0.5 * beta},{IMAG, 0.9 * beta}};

    double dt = 0.01;
    double dtau = 0.01;


    int N_t = int(t_max / dt + 0.1) + 1;
    int N_tau = int(beta / dtau + 0.1)  + 1;
    std::cout << N_tau  << std::endl;

    DysonEquation dyson_eq(t_max, beta, N_t, N_tau, eps_0, eps, U_0, U, V_0, V, eps_bath_0, eps_bath);
    dyson_eq.generate_du(times, n_ini);
    dyson_eq.generate_G0();
    dyson_eq.generate_F_and_M();
    dyson_eq.solve_for_G();

    std::cout << dyson_eq.Z() << ", " << dyson_eq.Z_from_G() << std::endl;

    arma::cx_mat G = dyson_eq.get_G() + arma::diagmat(dyson_eq.get_G_jump());

    int n_t_max = N_t - 1;
    std::cout << complex_t{0, 1} * G(n_t_max, n_t_max) << std::endl;


    return 0;








    std::vector<double> error;
    std::vector<complex_t> Z;
    std::vector<complex_t> Z_from_G;


    std::vector<double> ds_vec;
    for(int i = 51; i < 402; i += 10){
        ds_vec.push_back(beta / (i - 1));
    }
    for(int i = 56; i < 407; i += 10){
        ds_vec.push_back(beta / (i - 1));
    }
    for(int i = 53; i < 404; i += 10){
        ds_vec.push_back(beta / (i - 1));
    }
    for(auto ds : ds_vec){
        dt = ds;
        dtau = ds;
        int N_t = 0;// int(t_max / dt + 0.1) + 1;
        int N_tau = int(beta / dtau + 0.1)  + 1;
        std::cout << "(N_t, N_tau) = " << N_t << ", " << N_tau << std::endl;

        int N = 2 * N_t + N_tau - (N_t > 0 ? 1 : 0) - (N_t * N_tau > 0 ? 1 : 0);

        DysonEquation dyson_eq(dt, dtau, N_t, N_tau, eps_0, eps, U_0, U, V_0, V, eps_bath_0, eps_bath);
        for(int i = 0; i < 1; ++i) {
            dyson_eq.generate_du(times, n_ini);
            dyson_eq.generate_G0();
            dyson_eq.generate_F_and_M();
            dyson_eq.solve_for_G();
        }

        matrix_t G = dyson_eq.get_G();
        vector_t G_jump = dyson_eq.get_G_jump();

        matrix_t BC(2, 2);
        BC(0, 0) = G(0, 0) - G_jump(0);
        BC(1, 0) = -G(N - 1, 0);
        BC(0, 1) = -G(0, N - 1) - 2. * G_jump(N - 1);
        BC(1, 1) = G(N - 1, N - 1) - G_jump(N - 1);

        std::cout << "ds = " << ds << std::endl;
        std::cout << BC << std::endl;
        double errors = 0;
        for(int i = 0; i < 2; ++i){
            for(int j = 0; j < 2; ++j){
                for(int k = 0; k < 2; ++k){
                    for(int l = 0; l < 2; ++l){
                        errors += std::abs(BC(i ,j) - BC(k, l));
                    }
                }
            }
        }
        error.push_back(errors / 2. / 6.);

        std::cout << dyson_eq.Z() << " " << dyson_eq.Z_from_G() << std::endl;
        Z.push_back(dyson_eq.Z());
        Z_from_G.push_back(dyson_eq.Z_from_G());
    }

    Gnuplot gp;
    arma::mat data;

    data.zeros(ds_vec.size(), 2);
    data.col(0) = arma::conv_to<arma::colvec>::from(ds_vec);
    data.col(1) = arma::conv_to<arma::colvec>::from(error);

    gp << "set xrange [0:0.025]" << std::endl;
    gp << "set yrange [0:*]" << std::endl;
    gp << "plot '-' using 1:2 with points title 'error(ds)' " << std::endl;
    gp.send1d(data);

    std::cin.get();

    data.zeros(ds_vec.size(), 3);
    data.col(0) = arma::conv_to<arma::colvec>::from(ds_vec);
    data.col(1) = arma::abs(arma::conv_to<arma::cx_colvec>::from(Z));
    data.col(2) = arma::abs(arma::conv_to<arma::cx_colvec>::from(Z_from_G));


    gp << "set xrange [0:0.025]" << std::endl;
    gp << "set yrange [*:*]" << std::endl;
    gp << "plot '-' using 1:2 with points title 'Z', '-' using 1:3 with points title 'Z from G'" << std::endl;
    gp.send1d(data);
    gp.send1d(data);







//    if(N_t > 0){
//        data.zeros(N_t, 5);
//        for(int i = 0; i < N_t; ++i){
//            data.col(0)(i) = i * dt;
//            data.col(1)(i) = std::IMAG(G(i + 1, i));
//            data.col(2)(i) = std::IMAG(G(i, i + 1));
//            data.col(3)(i) = std::IMAG(G(i, i) + G_jump(i));
//            data.col(4)(i) = std::IMAG(G(i, i) - G_jump(i));
//        }
//
//        gp << "set yrange [-1.1:1.1]" << std::endl;
//        gp << "plot '-' using 1:2 with points title 'G(t+1, t)', '-' using 1:3 with points title 'G(t, t+1)',"
//           <<      "'-' using 1:4 with points title 'G(t+0, t)', '-' using 1:5 with points title 'G(t, t+0)'" << std::endl;
//        gp.send1d(data);
//        gp.send1d(data);
//        gp.send1d(data);
//        gp.send1d(data);
//
//        std::cin.get();
//    }
//
//    if(N_tau > 0){
//        data.zeros(N_tau, 2);
//        for(int i = 0; i < N_tau; ++i){
//            data.col(0)(i) = i * dtau;
//            data.col(1)(i) = std::IMAG(G(N - N_tau + i, N - N_tau) + (i == 0 ? G_jump(N - N_tau) : 0));
//        }
//
//        std::cout << 1 + data.col(1)(0) << " " <<  -data.col(1)(N_tau - 1) << std::endl;
//
//        gp << "set yrange [-1.1:0.1]" << std::endl;
//        gp << "plot '-' using 1:2 with points title 'G(tau+0, 0)'" << std::endl;
//        gp.send1d(data);
//    }

}
