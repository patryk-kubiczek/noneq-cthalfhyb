#include "DiscreteBath.h"

using namespace arma;

DiscreteBath::DiscreteBath(double eps_0, double eps, double U_0, double U, const std::vector<cx_double> &V_0,
                           const std::vector<cx_double> &V, const std::vector<double> &eps_bath_0,
                           const std::vector<double> &eps_bath, double beta, double t_max, int n_blocks)
        : eps_0(eps_0), eps(eps), U_0(U_0), U(U), beta(beta), t_max(t_max), n_blocks(n_blocks){
    generate_e_and_O(V_0, V, eps_bath_0, eps_bath);
}

void DiscreteBath::generate_e_and_O(const std::vector<cx_double> &V_0, const std::vector<cx_double> &V,
                                    const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath) {



    N = V.size() + 1;

    // R: real-time branch, I: imaginary-time branch

    // Imaginary time evolution

    mat E_I(N, N, fill::zeros);
    O_I_0.zeros(N, N);
    O_I_1.zeros(N, N);

    E_I(0, 0) = eps_0;

    for (int i = 1; i < N; ++i) {
        E_I(0, i) = std::abs(V_0[i - 1]);
        E_I(i, 0) = std::abs(V_0[i - 1]);
        E_I(i, i) = eps_bath_0[i - 1];
    }

    eig_sym(e_I_0, O_I_0, E_I);

    E_I(0, 0) += U_0;
    eig_sym(e_I_1, O_I_1, E_I);

    // Real time evolution

    cx_mat E(N, N, fill::zeros);
    O_R_0.zeros(N, N);
    O_R_1.zeros(N, N);

    E(0, 0) = eps;
    for (int i = 1; i < N; ++i) {
        E(0, i) = std::conj(V[i - 1]);
        E(i, 0) = V[i - 1];
        E(i, i) = eps_bath[i - 1];
    }

    eig_sym(e_R_0, O_R_0, E);
    E(0, 0) += U;


    eig_sym(e_R_1, O_R_1, E);



    // Initialize transition matrices

    O_I_01 = O_I_0.t() * O_I_1;
    O_R_01 = O_R_0.t() * O_R_1;

    O_IR_0 = O_I_0.t() * O_R_0;
    O_IR_1 = O_I_1.t() * O_R_1;

    // Initialize the matrices/vectors

    Q.zeros(N, N, n_blocks);
    R.zeros(N, N, n_blocks);
    D.zeros(N, n_blocks);
    T.zeros(N, N, n_blocks);

    D_b.zeros(N);
    D_s.zeros(N);

    u_t.zeros(N, N, 2);
    u_tau.zeros(N, N, n_blocks);

}

void DiscreteBath::generate_u(const std::vector<contour_time> &time_points, int n_ini) {

    const double dtau = beta / n_blocks;

    int n = n_ini;
    int i = 0;
    double prev_time = 0;
    double next_time = 0;
    // Forward branch
    u_t.slice(0) = (n == 0 ? O_IR_0.t() : O_IR_1.t());
    while (i < time_points.size() && time_points[i].branch == PLUS) {
        next_time = time_points[i].t;
        //DEBUG( std::cout << "Forward branch, prev_time: " << prev_time << ", next_time: " << next_time << ", n:" << n << std::endl; )
        if(n == 0)
            u_t.slice(0) = O_R_01.t() * diagmat(exp(-I * (next_time - prev_time) * e_R_0)) * u_t.slice(0);
        else
            u_t.slice(0) = O_R_01 * diagmat(exp(-I * (next_time - prev_time) * e_R_1)) * u_t.slice(0);
        prev_time = next_time;
        n = (n + 1) % 2;
        ++i;
    }
    next_time = t_max;
    //DEBUG( std::cout << "Forward branch, prev_time: " << prev_time << ", next_time: " << next_time << ", n:" << n << std::endl; )
    u_t.slice(0) = diagmat(exp(-I * (next_time - prev_time) * (n == 0 ? e_R_0 : e_R_1))) * u_t.slice(0);

    // Remember n at t_max
    n_at_t_max = n;

    //Backward branch
    u_t.slice(1).eye();
    prev_time = t_max;
    while (i < time_points.size() && time_points[i].branch == MINUS) {
        next_time = time_points[i].t;
        //DEBUG( std::cout << "Backward branch, prev_time: " << t_max - prev_time << ", next_time: " << t_max - next_time << ", n:" << n << std::endl;  )
        if(n == 0)
            u_t.slice(1) = O_R_01.t() * diagmat(exp(-I * (next_time - prev_time) * e_R_0)) * u_t.slice(1);
        else
            u_t.slice(1) = O_R_01 * diagmat(exp(-I * (next_time - prev_time) * e_R_1)) * u_t.slice(1);
        prev_time = next_time;
        n = (n + 1) % 2;
        ++i;
    }
    next_time = 0;
    //DEBUG( std::cout << "Backward branch, prev_time: " << t_max - prev_time << ", next_time: " << t_max - next_time << ", n:" << n << std::endl; )
    if(n == 0)
        u_t.slice(1) = O_IR_0 * diagmat(exp(-I * (next_time - prev_time) * e_R_0)) * u_t.slice(1);
    else
        u_t.slice(1) = O_IR_1 * diagmat(exp(-I * (next_time - prev_time) * e_R_1)) * u_t.slice(1);

    // Imaginary branch
    for (int j = 0; j < n_blocks; ++j) {
        u_tau.slice(j).eye();
        prev_time = j * dtau;
        while (i < time_points.size() && time_points[i].t < (j + 1) * dtau) {
            next_time = time_points[i].t;
            //DEBUG( std::cout << "Imaginary branch, block: " << j << ", prev_time: " <<  prev_time << ", next_time: " << next_time << ", n:" << n << std::endl; )
            if(n == 0)
                u_tau.slice(j) = O_I_01.t() * diagmat(exp(-(next_time - prev_time) * e_I_0)) * u_tau.slice(j);
            else
                u_tau.slice(j) = O_I_01 * diagmat(exp(-(next_time - prev_time) * e_I_1)) * u_tau.slice(j);
            prev_time = next_time;
            n = (n + 1) % 2;
            ++i;
        }
        next_time = (j + 1) * dtau;
        //DEBUG( std::cout << "Imaginary branch, block: " << j << ", prev_time: " <<  prev_time << ", next_time: " << next_time << ", n:" << n << std::endl; )
        u_tau.slice(j) = diagmat(exp(-(next_time - prev_time) * (n == 0 ? e_I_0 : e_I_1))) * u_tau.slice(j);
    }
}

void DiscreteBath::generate_QDT() {

    qr(Q.slice(0), R.slice(0), u_tau.slice(0));
    D.col(0) = abs(R.slice(0).diag());
    T.slice(0) = inv(diagmat(D.col(0))) * trimatu(R.slice(0));

    for(int i = 1; i < n_blocks; ++i){
        qr(Q.slice(i), R.slice(i), u_tau.slice(i) * Q.slice(i - 1) * diagmat(D.col(i - 1)));
        D.col(i) = abs(R.slice(i).diag());
        T.slice(i) = inv(diagmat(D.col(i))) * trimatu(R.slice(i)) * trimatu(T.slice(i - 1));
        //std::cout << D.col(i) << std::endl;
        //std::cout << T.slice(i) << std::endl;
    }

    D_b.ones();
    D_s.ones();
    for(int i = 0; i < N; ++i){
        if(D(i, n_blocks - 1) > 1)
            D_b(i) = D(i, n_blocks - 1);
        else
            D_s(i) = D(i, n_blocks - 1);
    }

}

cx_double DiscreteBath::log_z() {

    // determinant = det(1 + Q D T U^-_t U^+_t) = det(1 + Q D_b D_s T U^-_t U^+_t) = det(Q D_b) det ( D_b^-1 Q^-1 + D_s T U^-_t U^+_t)

    return  log_det(diagmat(D_b)) + log_det(Q.slice(n_blocks - 1))
            + log_det(inv(diagmat(D_b)) * Q.slice(n_blocks - 1).t()
                      + diagmat(D_s) * trimatu(T.slice(n_blocks - 1)) * u_t.slice(1) * u_t.slice(0));
}

cx_double DiscreteBath::log_z_without_QDT() {

    arma::mat u_tau_multiplied(u_tau.slice(0));
    for(int i = 1; i < n_blocks; ++i)
        u_tau_multiplied = u_tau.slice(i) * u_tau_multiplied;

    return log_det(diagmat(ones<cx_vec>(N)) + u_tau_multiplied * u_t.slice(1) * u_t.slice(0));
}

cx_mat DiscreteBath::density_matrix_at_t_max() const {

    // (1 + U^+_t Q D T U^-_t) G = 1 <-> (1 + U^+_t Q D_b D_s T U^-_t) G = 1
    // <-> (D_b^-1  Q^-1 U^+_t^-1  +  D_s T U^-_t) G = D_b^-1  Q^-1 U^+_t^-1  <-> A G = B

    cx_mat B = inv(diagmat(D_b)) * Q.slice(n_blocks - 1).t() * u_t.slice(0).t();
    cx_mat G = solve(B + diagmat(D_s) * trimatu(T.slice(n_blocks - 1)) * u_t.slice(1), B);

    return (n_at_t_max == 0 ? O_R_0 * G * O_R_0.t() : O_R_1 * G * O_R_1.t());
}


cx_mat DiscreteBath::density_matrix_at_t_max_without_QDT() const {

    arma::mat u_tau_multiplied(u_tau.slice(0));
    for(int i = 1; i < n_blocks; ++i)
        u_tau_multiplied = u_tau.slice(i) * u_tau_multiplied;

    cx_mat G = inv(diagmat(ones<cx_vec>(N)) + u_t.slice(0) * u_tau_multiplied * u_t.slice(1));

    return (n_at_t_max == 0 ? O_R_0 * G * O_R_0.t() : O_R_1 * G * O_R_1.t());
}




