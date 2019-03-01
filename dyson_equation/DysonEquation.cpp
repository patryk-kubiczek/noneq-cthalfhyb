#include "DysonEquation.h"

using namespace arma;

DysonEquation::DysonEquation(double t_max, double beta, int N_t, int N_tau,
                             double eps_0, double eps, double U_0, double U,
                             const cx_mat &Delta,  const cx_vec &Delta_jump)
        :
        dt(t_max / (N_t - 1)), dtau(beta / (N_tau - 1)), N_t(N_t), N_tau(N_tau),
        eps_0(eps_0), eps(eps), U_0(U_0), U(U),
        Delta(Delta), Delta_jump(Delta_jump)  { initialize(); }

DysonEquation::DysonEquation(double t_max, double beta, int N_t, int N_tau, double eps_0, double eps, double U_0, double U,
                             const std::vector<cx_double> &V_0, const std::vector<cx_double> &V,
                             const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath)
        :
        dt(t_max / (N_t - 1)), dtau(beta / (N_tau - 1)), N_t(N_t), N_tau(N_tau),
        eps_0(eps_0), eps(eps), U_0(U_0), U(U) {
    initialize();
    assert(V.size() == V_0.size() && V.size() == eps_bath.size() && V.size() == eps_bath_0.size());
    generate_simple_Delta(V_0, V, eps_bath_0, eps_bath);
}


void DysonEquation::initialize() {
    assert(N_t == 0 || N_t > 1);
    assert(N_tau == 0 || N_tau > 1);
    assert(N_t > 0 || N_tau > 0);
    N = 2 * N_t + N_tau;
    if(N_t > 0){
        --N;
        i_t_max = N_t - 1;
        i_0 = 2 * N_t - 2;
    }
    else{
        i_t_max = 0;
        i_0 = 0;
    }
    if(N_t > 0 && N_tau > 0) --N;

    du.resize(N - 1);
    G0.zeros(N, N);
    G0_jump.zeros(N);
    F.zeros(N, N);
    M.zeros(N, N);
    G.zeros(N, N);
    generate_w();

}

void DysonEquation::generate_w() {
    W.zeros(N);
    if(N_t > 0) W(0) = dt / 2.;
    for(int i = 1; i < i_t_max; ++i)
        W(i) = dt;
    if(N_t > 0) W(i_t_max) = 0;
    for(int i = i_t_max + 1; i < i_0; ++i)
        W(i) = -dt;
    if(N_t > 0) W(i_0) = -dt / 2.;
    if(N_tau > 0) W(i_0) += cx_double{0, -1} * dtau / 2.;
    for(int i = i_0 + 1; i < N - 1; ++i)
        W(i) = cx_double{0, -1} * dtau;
    if(N_tau > 0) W(N - 1) = cx_double{0, -1} * dtau / 2.;
}



void DysonEquation::generate_du(const std::vector<contour_time>& time_points, int n_ini) {
    cx_double du0_t = std::exp(-I * dt * eps);
    cx_double du0_tau = std::exp(-dtau * eps_0);

    cx_double du1_t = std::exp(-I * dt * (eps + U));
    cx_double du1_tau = std::exp(-dtau * (eps_0 + U_0));

    auto dux_t = [dt=dt, eps=eps, U=U](double x) -> cx_double {
        return std::exp(-I * dt * (eps + x * U));
    };
    auto dux_tau = [dtau=dtau, eps_0=eps_0, U_0=U_0](double x) -> cx_double {
        return std::exp(-dtau * (eps_0 + x * U_0));
    };

    if(time_points.empty()){
        cx_double du_t = (n_ini == 0 ? du0_t : du1_t);
        cx_double du_tau = (n_ini == 0 ? du0_tau : du1_tau);
        for(int i = 0; i < i_t_max; ++i) du[i]  = du_t;
        du_t = std::conj(du_t);
        for(int i = i_t_max; i < i_0; ++i) du[i] = du_t;
        for(int i = i_0; i < N - 1; ++i) du[i] = du_tau;
        n_down_at_t_max = n_ini;
    }
    else{
        int n = n_ini;
        int i_change = 0;
        contour_time t_change = time_points[i_change];
        cx_double du_t = (n_ini == 0 ? du0_t : du1_t);
        for(int i = 0; i < i_t_max; ++i){
            if(t_change.branch == PLUS && (i + 1) * dt > t_change.t){
                double x = (n == 0 ? i + 1 - t_change.t / dt : t_change.t / dt - i);
                du[i] = dux_t(x);
                n = (n + 1) % 2;
                t_change = (++i_change < time_points.size() ? time_points[i_change] : contour_time{IMAG, N_tau * dtau});
                du_t = (n == 0 ? du0_t : du1_t);
            }
            else du[i] = du_t;
        }
        n_down_at_t_max = n;
        du_t = std::conj(du_t);
        for(int i = i_t_max; i < i_0; ++i){
            if(t_change.branch == MINUS && (i_t_max - (i - i_t_max + 1)) * dt < t_change.t){
                double x = (n == 0 ? -(i_t_max - (i - i_t_max + 1)) + t_change.t / dt
                                   : i_t_max - (i - i_t_max) - t_change.t / dt );
                du[i] = std::conj(dux_t(x));
                n = (n + 1) % 2;
                t_change = (++i_change < time_points.size() ? time_points[i_change] : contour_time{IMAG, N_tau * dtau});
                du_t = std::conj((n == 0 ? du0_t : du1_t));
            }
            else du[i] = du_t;
        }
        cx_double du_tau = (n == 0 ? du0_tau : du1_tau);
        for(int i = i_0; i < N - 1; ++i){
            if(t_change.branch == IMAG && (i - i_0 + 1) * dtau > t_change.t){
                double x = (n == 0 ? i - i_0 + 1 - t_change.t / dtau : t_change.t / dtau - i + i_0);
                du[i] = dux_tau(x);
                n = (n + 1) % 2;
                t_change = (++i_change < time_points.size() ? time_points[i_change] : contour_time{IMAG, N_tau * dtau});
                du_tau = (n == 0 ? du0_tau : du1_tau);
            }
            else du[i] = du_tau;
        }
    }
}

void DysonEquation::generate_simple_Delta(const std::vector<cx_double> &V_0, const std::vector<cx_double> &V,
                                          const std::vector<double> &eps_bath_0, const std::vector<double> &eps_bath) {
    const int N_bath = V.size();
    Delta.zeros(N, N);
    Delta_jump.zeros(N);
    std::vector<cx_double> du(N - 1);
    cx_mat g(N, N);

    for(int n = 0; n < N_bath; ++n){

        cx_double du_t = std::exp(-I * dt * eps_bath[n]);
        cx_double du_tau = std::exp(-dtau * eps_bath_0[n]);
        for(int i = 0; i < i_t_max; ++i)
            du[i] = du_t;
        for(int i = i_t_max; i < i_0; ++i)
            du[i] = std::conj(du_t);
        for(int i = i_0; i < N - 1; ++i)
            du[i] = du_tau;

        for(int j = 0; j < N; ++j){
            g(j, j) = 1.;
            for(int i = j + 1; i < N; ++i){
                g(i, j) = g(i - 1, j) * du[i - 1];
            }
            for(int i = j - 1; i >= 0; --i){
                g(i, j) = g(i + 1, j) / du[i];
            }
        }

        cx_double greater_prefix = -I / (1. + g(N - 1, 0));
        cx_double lesser_prefix = I / (1. + g(0, N - 1));
        cx_double mean = (greater_prefix + lesser_prefix) / 2.;
        cx_double difference = (greater_prefix - lesser_prefix) / 2.;

        Z_bath *= 1. + g(N - 1, 0);

        for(int j = 0; j < N; ++j){
            for(int i = 0; i < j; ++i)
                g(i, j) *= lesser_prefix;
            g(j, j) = mean;
            for(int i = j + 1; i < N; ++i)
                g(i, j) *= greater_prefix;
        }
        for(int j = 0; j < i_0; ++j){
            for(int i = 0; i < i_0; ++i)
                Delta(i, j) += std::conj(V[n]) * g(i, j) * V[n];
            for(int i = i_0; i < N; ++i)
                Delta(i, j) += std::conj(V_0[n]) * g(i, j) * V[n];
        }
        for(int j = i_0; j < N - 1; ++j){
            for(int i = 0; i < i_0; ++i)
                Delta(i, j) += std::conj(V[n]) * g(i, j) * V_0[n];
            for(int i = i_0; i < N; ++i)
                Delta(i, j) += std::conj(V_0[n]) * g(i, j) * V_0[n];
        }
        for(int i = 0; i < i_0; ++i)
            Delta_jump(i) += V[n] * difference * V[n];
        for(int i = i_0; i < N; ++i)
            Delta_jump(i) += V_0[n] * difference * V_0[n];

    }



}

void DysonEquation::generate_G0() {
    for(int j = 0; j < N; ++j){
        G0(j, j) = 1.;
        for(int i = j + 1; i < N; ++i){
            G0(i, j) = G0(i - 1, j) * du[i - 1];
        }
        for(int i = j - 1; i >= 0; --i){
            G0(i, j) = G0(i + 1, j) / du[i];
        }
    }
    cx_double greater_prefix = -I / (1. + G0(N - 1, 0));
    cx_double lesser_prefix = I / (1. + G0(0, N - 1));
    cx_double mean = (greater_prefix + lesser_prefix) / 2.;
    cx_double difference = (greater_prefix - lesser_prefix) / 2.;

    Z0 = 1. + G0(N - 1, 0);

    for(int j = 0; j < N; ++j){
        for(int i = 0; i < j; ++i)
            G0(i, j) *= lesser_prefix;
        G0(j, j) = mean;
        for(int i = j + 1; i < N; ++i)
            G0(i, j) *= greater_prefix;
    }
    for(int i = 0; i < N; ++i)
        G0_jump(i) = difference;


}

matrix_t DysonEquation::convolution(const cx_mat& A, const cx_mat &B, const cx_vec& DA, const cx_vec& DB){
    cx_mat C = A * diagmat(W) * B;

    C.diag() -= W % DA % DB;
    C.row(0) -= W(0) * DA(0) * B.row(0);
    C.col(0) += W(0) * A.col(0) * DB(0);
    C.row(N - 1) += W(N - 1) * DA(N - 1) * B.row(N - 1);
    C.col(N - 1) -= W(N - 1) * A.col(N - 1) * DB(N - 1);

    if(N_t > 0){
        C.row(i_t_max) += dt * DA(i_t_max) * B.row(i_t_max);
        C.col(i_t_max) -= dt * A.col(i_t_max) * DB(i_t_max);
    }
    if(N_t > 0 && N_tau > 0) {
        C.row(i_0) += (-dt + I * dtau) / 2. * DA(i_0) * B.row(i_0);
        C.col(i_0) -= (-dt + I * dtau) / 2. * A.col(i_0) * DB(i_0);
    }
    return C;
}

void DysonEquation::generate_F_and_M() {
    F = -convolution(G0, Delta, G0_jump, Delta_jump);
    M = F * diagmat(W);
    M.diag() += 1.;
}


void DysonEquation::solve_for_G() {
    // (1 - G0 o Delta) o G = (1 - F) o G = M * G - \delta(F o G) = G0
    cx_mat M_inv = inv(M);
    G = solve(M, G0, solve_opts::fast);
    cx_mat delta_G(N, N);
    int i = 0;
    while(i < 2){
        delta_G = M_inv * (G0 - G - convolution(F, G, cx_vec(N, fill::zeros), G0_jump));
        std::cout << "Norm(delta) = " << norm(delta_G) << std::endl;
        G += delta_G;
        ++i;
    }
}

cx_double DysonEquation::Z() {
    //Z0 is calculated in generate_G0()
    return Z_bath * Z0 * det(M);
}

cx_double DysonEquation::Z_from_G() {
    //Z0 is calculated in generate_G0()
    return Z_bath / det(I * (G + diagmat(G0_jump)));
}











