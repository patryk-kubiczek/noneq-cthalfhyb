#include "InputGenerator.h"
#include "../auxiliary_functions/constants.h"

using namespace arma;

InputGenerator::InputGenerator(const UserInput& data) : data(data){
    // Make sure hybridizations and bath-energies have the same size
    assert(data.bath_energies.size() == data.hybridizations.size());

    // Make sure the number of time points is the same in each data vector
    for(const auto &k_vec : data.bath_energies) {
        for(const auto &flavor_vec : k_vec) {
            assert(flavor_vec.size() == data.n_times);
        }
    }
    for(const auto &k_vec : data.hybridizations) {
        for(const auto &flavor_vec : k_vec) {
            assert(flavor_vec.size() == data.n_times);
        }
    }
    for(const auto &flavor_vec : data.local_energies) {
        assert(flavor_vec.size() == data.n_times);
    }
    assert(data.Hubbard_U.size() == data.n_times);

    // Make sure the number of flavors is the same
    int n_flavor = data.local_energies.size();
    assert(n_flavor == data.bath_energies[0].size() && n_flavor == data.hybridizations[0].size());
}

void InputGenerator::generate_p_t_t() {

    std::cout << "Generating atomic propagators on real axis" << std::endl;

    // This implementation is valid only for n_flavor <= 2
    int n_flavor = data.local_energies.size();
    assert(n_flavor <= 2);

    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate H_loc(t)
    cx_cube H_loc;
    if(n_flavor == 1) {
        const vec_re_1_t &energy = data.local_energies[0];
        H_loc.zeros(2, 2, data.n_times);
        H_loc.tube(1, 1) = conv_to<cx_vec>::from(energy);
    }
    if(n_flavor == 2) {
        const vec_re_1_t &energy_up = data.local_energies[0];
        const vec_re_1_t &energy_down = data.local_energies[1];
        H_loc.zeros(4, 4, data.n_times);
        H_loc.tube(1, 1) = conv_to<cx_vec>::from(energy_up);
        H_loc.tube(2, 2) = conv_to<cx_vec>::from(energy_down);
        H_loc.tube(3, 3) = conv_to<cx_vec>::from(energy_up)
                           + conv_to<cx_vec>::from(energy_down)
                           + conv_to<cx_vec>::from(data.Hubbard_U);
    }

    EquationsOfMotion<> eom(H_loc, data.dt, data.t_max, data.real_grid_size);
    eom.run();
    double error = eom.check_unitarity();
    std::cout << "Testing unitarity, maximal error: " << error << std::endl;

    //double real_dt = data.t_max / (data.real_grid_size - 1);
    eom.save(get_filename(data.params_name, p_t_t_name));
}


void InputGenerator::generate_p_tau_tau() {

    std::cout << "Generating atomic propagators on imaginary axis" << std::endl;

    int n_flavor = data.local_energies.size();
    assert(n_flavor <= 2);

    std::cout << "Number of flavors: " << n_flavor << std::endl;

    // Generate H_loc(0)
    mat H_loc_0;
    if(n_flavor == 1) {
        const vec_re_1_t &energy = data.local_energies[0];
        H_loc_0.zeros(2, 2);
        H_loc_0(1, 1) = energy[0];
    }
    if(n_flavor == 2) {
        const vec_re_1_t &energy_up = data.local_energies[0];
        const vec_re_1_t &energy_down = data.local_energies[1];
        H_loc_0.zeros(4, 4);
        H_loc_0(1, 1) = energy_up[0];
        H_loc_0(2, 2) = energy_down[0];
        H_loc_0(3, 3) = energy_up[0] + energy_down[0] + data.Hubbard_U[0];
    }

    cube p(H_loc_0.n_rows, H_loc_0.n_cols, data.imag_grid_size);

    double imag_dt = data.beta / (data.imag_grid_size - 1);

    for(int i = 0; i < p.n_slices; ++i) {
        p.slice(i) = expmat(-i * imag_dt * H_loc_0);
    }

    p.save(get_filename(data.params_name, p_tau_tau_name));
}


void InputGenerator::generate_delta() {

    std::cout << "Generating hybridization functions" << std::endl;

    int n_k = data.bath_energies.size();
    int n_flavor = data.bath_energies[0].size();

    std::cout << "Number of bath sites: " << n_k << std::endl;
    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Finding evolution operator on the real branch for each k ***********************************
    field<cx_tetracube> U_bath_t_t_k(n_k); // output of eom goes here
    U_bath_t_t_k.fill(cx_tetracube(n_flavor, n_flavor, data.real_grid_size, data.real_grid_size));
    cx_cube H_bath(n_flavor, n_flavor, data.n_times, fill::zeros);

    double max_error = 0;
    double error = 0;

    // Assuming no off-diagonal (in flavor space) H_bath elements
    for(int k = 0; k < n_k; ++k) {
        H_bath.zeros();
        for(int a = 0; a < n_flavor; ++a) {
            H_bath.tube(a, a) = conv_to<cx_vec>::from(data.bath_energies[k][a]); // Setting diagonal elements
        }
        EquationsOfMotion<> eom(H_bath, data.dt, data.t_max, data.real_grid_size);
        eom.run();
        error = eom.check_unitarity();
        if(error > max_error){
            max_error = error;
        }
        U_bath_t_t_k(k) = eom.get_result();
    }
    std::cout << "Testing unitarity, maximal error: " << max_error << std::endl;
    // ********************************************************************************************


    // Finding evolution operator on the imaginary branch for each k ******************************
    field<cx_cube> U_bath_tau_tau_k(n_k); // output of exponentiation goes here
    U_bath_tau_tau_k.fill(cx_cube(n_flavor, n_flavor, data.imag_grid_size, fill::zeros));

    cx_mat H_bath_0(n_flavor, n_flavor);
    double imag_dt = (data.imag_grid_size > 0 ? data.beta / (data.imag_grid_size - 1) : 0);

    // Assuming no off-diagonal (in flavor space) H_bath elements
    for(int k = 0; k < n_k; ++k) {
        H_bath_0.zeros();
        for(int a = 0; a < n_flavor; ++a) {
            H_bath_0(a, a) = data.bath_energies[k][a][0]; // Setting diagonal elements
        }
        cx_cube U(n_flavor, n_flavor, data.imag_grid_size);
        for(int i = 0; i < U.n_slices; ++i) {
            U.slice(i) = expmat(-i * imag_dt * diagmat(H_bath_0));
        }
        U_bath_tau_tau_k(k) = U;
    }
    // ********************************************************************************************



    // Creating greater (>) and lesser (<) t-t Green function ********************************************************
    field<cx_tetracube> g_bath_greater_k(n_k);
    g_bath_greater_k.fill(cx_tetracube(n_flavor, n_flavor, data.real_grid_size, data.real_grid_size));

    field<cx_tetracube> g_bath_lesser_k(n_k);
    g_bath_lesser_k.fill(cx_tetracube(n_flavor, n_flavor, data.real_grid_size, data.real_grid_size));

    cx_mat hole_fermi_function(n_flavor, n_flavor, fill::zeros);
    cx_mat fermi_function(n_flavor, n_flavor, fill::zeros);

    for(int k = 0; k < n_k; ++k) {
        cx_mat& U_beta_0  = U_bath_tau_tau_k(k).slice(data.imag_grid_size - 1);
        hole_fermi_function = inv(diagmat(eye<cx_mat>(n_flavor, n_flavor) + U_beta_0));
        fermi_function = diagmat(hole_fermi_function) * U_beta_0;

        for(int i = 0; i < data.real_grid_size; ++i){
            for(int j = 0; j < data.real_grid_size; ++j){
                g_bath_greater_k(k).slice(i, j) =  -I * U_bath_t_t_k(k).slice(i, 0) * hole_fermi_function * U_bath_t_t_k(k).slice(0, j);
                g_bath_lesser_k(k).slice(i, j) =  I * U_bath_t_t_k(k).slice(i, 0) * fermi_function * U_bath_t_t_k(k).slice(0, j);
            }
        }
    }
    // ********************************************************************************************

    // Creating tau-t and tau-tau Green function **************************************************************
    field<cx_tetracube> g_bath_tau_t_k(n_k);
    g_bath_tau_t_k.fill(cx_tetracube(n_flavor, n_flavor, data.imag_grid_size, data.real_grid_size));

    // g_bath_tau_tau not needed since g_bath(tau) = g_bath(tau, t=0)
    //field<matrix_storage<1>> g_bath_tau_tau_k(n_k);
    //g_bath_tau_tau_k.fill(matrix_storage<1>(n_flavor, n_flavor, data.imag_grid_size));

    for(int k = 0; k < n_k; ++k) {
        cx_mat& U_beta_0  = U_bath_tau_tau_k(k).slice(data.imag_grid_size - 1);
        hole_fermi_function = inv(diagmat(eye<cx_mat>(n_flavor, n_flavor) + U_beta_0));
        for(int i = 0; i < data.imag_grid_size; ++i){
            //g_bath_tau_tau_k(k)(i) =  -I * U_bath_tau_tau_k(k)(i) * hole_fermi_function;
            for(int j = 0; j < data.real_grid_size; ++j){
                g_bath_tau_t_k(k).slice(i, j) =  -I * diagmat(U_bath_tau_tau_k(k).slice(i))
                        * diagmat(hole_fermi_function) * diagmat(U_bath_t_t_k(k).slice(0, j));
            }
        }
    }
    // ********************************************************************************************


    // Generate deltas *************************************************************
    field<cx_mat> delta_greater_t_t(n_flavor, n_flavor);
    delta_greater_t_t.fill(cx_mat(data.real_grid_size, data.real_grid_size, fill::zeros));

    field<cx_mat> delta_lesser_t_t(n_flavor, n_flavor);
    delta_lesser_t_t.fill(cx_mat(data.real_grid_size, data.real_grid_size, fill::zeros));

    field<cx_mat> delta_tau_t(n_flavor, n_flavor);
    delta_tau_t.fill(cx_mat(data.imag_grid_size, data.real_grid_size, fill::zeros));

    //field<vector_t> delta_tau_tau(n_flavor, n_flavor);
    //delta_tau_tau.fill(vector_t(data.imag_grid_size));

    // Works only for diagonal hybridization
    for(int k = 0; k < n_k; ++k){
        for(int a = 0; a < n_flavor; ++a){
            for(int b = 0; b < n_flavor; ++b){
                for(int i = 0; i < data.real_grid_size; ++i){
                    for(int j = 0; j < data.real_grid_size; ++j){
                        auto &hyb_a = data.hybridizations[k][a][i];
                        auto &hyb_b = data.hybridizations[k][b][j];
                        delta_greater_t_t(a, b)(i, j) += std::conj(hyb_a) * g_bath_greater_k(k)(a, b, i, j) * hyb_b;
                        delta_lesser_t_t(a, b)(i, j) += std::conj(hyb_a) * g_bath_lesser_k(k)(a, b, i, j) * hyb_b;
                    }
                }
            }
        }
    }

    for(int k = 0; k < n_k; ++k){
        for(int a = 0; a < n_flavor; ++a){
            for(int b = 0; b < n_flavor; ++b){
                for(int i = 0; i < data.imag_grid_size; ++i){
                    auto &hyb_a = data.hybridizations[k][a][0];
                    //cx_double &hyb_b = data.hybridizations(k)(b)(0);
                    //delta_tau_tau(a, b)(i) = std::conj(hyb_a) * g_bath_tau_tau_k(k)(i)(a, b) * hyb_b;
                    for(int j = 0; j < data.real_grid_size; ++j){
                        auto &hyb_b = data.hybridizations[k][b][j];
                        delta_tau_t(a, b)(i, j) += std::conj(hyb_a) * g_bath_tau_t_k(k)(a, b, i, j) * hyb_b;
                    }
                }
            }
        }
    }

    delta_lesser_t_t.save(get_filename(data.params_name, delta_lesser_t_t_name));
    delta_greater_t_t.save(get_filename(data.params_name, delta_greater_t_t_name));
    delta_tau_t.save(get_filename(data.params_name, delta_tau_t_name));
    //delta_tau_tau.Saver(get_filename(data.params_name, delta_tau_tau_name));

    delta_lesser_t_t(0, 0).save(get_filename(data.params_name, delta_up_lesser_t_t_name));
    delta_greater_t_t(0, 0).save(get_filename(data.params_name, delta_up_greater_t_t_name));
    delta_tau_t(0, 0).save(get_filename(data.params_name, delta_up_tau_t_name));

    delta_lesser_t_t(1, 1).save(get_filename(data.params_name, delta_down_lesser_t_t_name));
    delta_greater_t_t(1, 1).save(get_filename(data.params_name, delta_down_greater_t_t_name));
    delta_tau_t(1, 1).save(get_filename(data.params_name, delta_down_tau_t_name));

//    delta_lesser_t_t(1, 1).save(get_filename(data.params_name, delta_down_lesser_t_t_name), arma_ascii);
//    delta_greater_t_t(1, 1).save(get_filename(data.params_name, delta_down_greater_t_t_name), arma_ascii);
//    delta_tau_t(1, 1).save(get_filename(data.params_name, delta_down_tau_t_name),  arma_ascii);

}

void InputGenerator::generate_u_t_t(int n) {
    std::cout << "Generating u_" << n << " on real axis" << std::endl;

    int spin_up_a = 0;
    // This input is only meaningful for n_flavor == 2
    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);

    int n_k = data.bath_energies.size();

    std::cout << "Number of bath sites: " << n_k << std::endl;
    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate h_up(t)
    cx_cube h_up(n_k + 1, n_k + 1, data.n_times, fill::zeros);
    h_up.tube(0, 0) = conv_to<cx_vec>::from(data.local_energies[spin_up_a]);

    if(n == 1){
        h_up.tube(0, 0) += conv_to<cx_vec>::from(data.Hubbard_U);

    }
    for(int k = 0; k < n_k; ++k){
        h_up.tube(k + 1, k + 1) = conv_to<cx_vec>::from(data.bath_energies[k][spin_up_a]);
        h_up.tube(0, k + 1) = conj(conv_to<cx_vec>::from(data.hybridizations[k][spin_up_a]));
        h_up.tube(k + 1, 0) = conv_to<cx_vec>::from(data.hybridizations[k][spin_up_a]);
    }


    EquationsOfMotion<> eom(h_up, data.dt, data.t_max, data.real_grid_size);
    eom.run();
    double error = eom.check_unitarity();
    std::cout << "Testing unitarity, maximal error: " << error << std::endl;

    //double real_dt = data.t_max / (data.real_grid_size - 1);
    eom.save(get_filename(data.params_name, (n == 0 ? u_0_t_t_name : u_1_t_t_name)));

}


void InputGenerator::generate_u_tau_tau(int n) {

    std::cout << "Generating u_" << n << " on imaginary axis" << std::endl;

    int spin_up_a = 0;
    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);
    int n_k = data.bath_energies.size();

    std::cout << "Number of bath sites: " << n_k << std::endl;
    std::cout << "Number of flavors: " << n_flavor << std::endl;

    // Generate h_up(0)
    mat h_up(n_k + 1, n_k + 1,  fill::zeros);
    h_up(0, 0) = data.local_energies[spin_up_a][0];
    if(n == 1)
        h_up(0, 0) += data.Hubbard_U[0];
    for(int k = 0; k < n_k; ++k){
        h_up(k + 1, k + 1) = data.bath_energies[k][spin_up_a][0];
        h_up(0, k + 1) = std::abs(std::conj(data.hybridizations[k][spin_up_a][0]));
        h_up(k + 1, 0) = std::abs(data.hybridizations[k][spin_up_a][0]);
    }

    cube u(h_up.n_rows, h_up.n_cols, data.imag_grid_size);

    if(data.imag_grid_size > 1) {

        double imag_dt = data.beta / (data.imag_grid_size - 1);

        for (int i = 0; i < u.n_slices; ++i) {
            mat m = -i * imag_dt * h_up;
            u.slice(i) = expmat_sym(m);

            mat X = eye<mat>(size(u.slice(i))) + u.slice(i);

            if (i == u.n_slices - 1) {
                //See whether this matrix is well-behaved
                std::cout << std::endl << "tau = " << i * imag_dt << std::endl;
                //std::cout << "Det from trace of exponent: " << std::exp(trace(m)) << std::endl;
                std::cout << "Det from det:               " << det(X) << std::endl;
                double val, sign;
                log_det(val, sign, X);
                std::cout << "Det from log_det:           " << std::exp(val) * sign << std::endl;

                mat U(size(X));
                mat V(size(X));
                vec s(X.n_cols);

                svd(U, s, V, u.slice(i), "std");

                mat Up(size(X));
                mat Vp(size(X));
                vec sp(X.n_cols);

                svd(Up, sp, Vp, inv(V.t() * U) + diagmat(s), "std");


                std::cout << "Det from SVD:               "
                          << det(U * Up) * det(diagmat(sp)) * det(Vp.t() * V.t()) << std::endl;

                vec eigvals;
                mat eigvecs;
                eig_sym(eigvals, eigvecs, X);

                double largest_eigenvalue_x = max(abs(eigvals));

                mat inverse_x = pinv(X);

                eig_sym(eigvals, eigvecs, inverse_x);

                double largest_eigenvalue_inv_x = max(abs(eigvals));

                std::cout << "Condition number: " << largest_eigenvalue_x * largest_eigenvalue_inv_x << std::endl
                          << std::endl;
            }
        }
    }
    else{
        u.slice(0) = expmat_sym(-data.beta * h_up);
    }


    u.save(get_filename(data.params_name, (n == 0 ? u_0_tau_tau_name : u_1_tau_tau_name)));
}


void InputGenerator::generate_phi_1_t_t() {
    std::cout << "Generating phi on real axis" << std::endl;

    int spin_down_a = 1;
    // This input is only meaningful for n_flavor == 2
    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);


    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate eps_down(t)
    cx_cube eps_down(1, 1, data.n_times,  fill::zeros);
    eps_down.tube(0, 0) = conv_to<cx_vec>::from(data.local_energies[spin_down_a]);

    EquationsOfMotion<> eom(eps_down, data.dt, data.t_max, data.real_grid_size);
    eom.run();
    double error = eom.check_unitarity();
    std::cout << "Testing unitarity, maximal error: " << error << std::endl;

    const cx_tetracube &result = eom.get_result();
    cx_mat phi(data.real_grid_size, data.real_grid_size, fill::zeros);

    for(int i = 0; i < result.n_slices_1(); ++i){
        for(int j = 0; j < result.n_slices_2(); ++j)
            phi(i, j) = result(0, 0, i, j);
    }

    //double real_dt = data.t_max / (data.real_grid_size - 1);
    phi.save(get_filename(data.params_name, phi_1_t_t_name));

}

void InputGenerator::generate_phi_1_tau_tau() {
    std::cout << "Generating phi on imaginary axis" << std::endl;

    int spin_down_a = 1;
    // This input is only meaningful for n_flavor == 2
    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);


    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate eps_down(0)
    double eps_down = data.local_energies[spin_down_a][0];

    vec phi(data.imag_grid_size);

    if(data.imag_grid_size > 1){
        double imag_dt = data.beta / (data.imag_grid_size - 1);

        for(int i = 0; i < phi.n_rows; ++i) {
            phi(i) = std::exp(-i * imag_dt * eps_down);
        }
    }
    else {
        phi(0) = std::exp(-data.beta * eps_down);
    }

    phi.save(get_filename(data.params_name, phi_1_tau_tau_name));
}

void InputGenerator::generate_c_and_a_operator_matrices() {
    std::cout << "Generating c and a operator matrices" << std::endl;

    // This is only valid for n_flavor == 2
    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);

    if(n_flavor == 2){
        cx_mat c_up = { {0, 0, 0, 0},
                          {1, 0, 0, 0},
                          {0, 0, 0, 0},
                          {0, 0, 1, 0} };
        cx_mat c_down = { {0, 0, 0, 0},
                            {0, 0, 0, 0},
                            {1, 0, 0, 0},
                            {0, -1, 0, 0} };
        cx_mat a_up = c_up.t();
        cx_mat a_down = c_down.t();

        cx_cube c_matrices(4, 4, 2);
        c_matrices.slice(0) = c_up;
        c_matrices.slice(1) = c_down;

        cx_cube a_matrices(4, 4, 2);
        a_matrices.slice(0) = a_up;
        a_matrices.slice(1) = a_down;

        c_matrices.save(get_filename(data.params_name, c_operator_matrices_name));
        a_matrices.save(get_filename(data.params_name, a_operator_matrices_name));
    }
}

void InputGenerator::generate_noninteracting_benchmark_data() {
    std::cout << "Generating benchmark data (U = 0, any n_bath_up)" << std::endl;

    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);

    int n_k = data.bath_energies.size();

    for(const auto& U : data.Hubbard_U){
        if(U != 0){
            std::cout << "Input U(t) != 0, exiting" << std::endl;
            return;
        }
    }


    std::cout << "Number of bath sites: " << n_k << std::endl;
    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate h(t)
    int N = n_flavor * (1 + n_k);
    int i = 0;
    cx_cube h(N, N, data.n_times, fill::zeros);
    for(int flavor = 0; flavor < n_flavor; ++flavor){
        i = flavor * (n_k + 1);
        h.tube(i, i) = conv_to<cx_vec>::from(data.local_energies[flavor]);
        for(int k = 0; k < n_k; ++k){
            h.tube(i + 1 + k, i + 1 + k) = conv_to<cx_vec>::from(data.bath_energies[k][flavor]);
            h.tube(i, i + 1 + k) = conj(conv_to<cx_vec>::from(data.hybridizations[k][flavor]));
            h.tube(i + 1 + k, i) = conv_to<cx_vec>::from(data.hybridizations[k][flavor]);
        }
    }

    std::cout << "h_t: " << std::endl << h.slice(data.n_times - 1) << std::endl;


    EquationsOfMotion<> eom(h, data.dt, data.t_max, data.real_grid_size);
    eom.run();
    double error = eom.check_unitarity();
    std::cout << "Testing unitarity, maximal error: " << error << std::endl;

    auto U_t_t = eom.get_result();


    cx_mat h_0(N, N, fill::zeros);

    for(int flavor = 0; flavor < n_flavor; ++flavor){
        i = flavor * (n_k + 1);
        h_0(i, i) = data.local_energies[flavor][0];
        for(int k = 0; k < n_k; ++k){
            h_0(i + 1 + k, i + 1 + k) = data.bath_energies[k][flavor][0];
            h_0(i, i + 1 + k) = std::conj(data.hybridizations[k][flavor][0]);
            h_0(i + 1 + k, i) = data.hybridizations[k][flavor][0];
        }
    }

    cx_mat fermi_function = inv(expmat_sym(data.beta * h_0) + eye<cx_mat>(size(h_0)));
    cx_mat hole_fermi_function = inv(expmat_sym(-data.beta * h_0) + eye<cx_mat>(size(h_0)));
    std::cout << "h_0: " << std::endl << h_0 << std::endl;
    std::cout << "fermi: " << std::endl << real(fermi_function) << std::endl;

    cx_cube U_tau(N, N, data.imag_grid_size);
    double imag_dt = data.beta / (data.imag_grid_size - 1);
    for(int j = 0; j < data.imag_grid_size; ++j) {
        U_tau.slice(j) = expmat_sym(-j * imag_dt * h_0);
    }

    int coarsing_factor = 1;
    int n_times = data.real_grid_size / coarsing_factor;
    int n_imag_times = data.imag_grid_size / coarsing_factor;

    std::cout << "imaginary time evolution matrix: " << std::endl << U_tau(n_imag_times - 1) << std::endl;

    std::cout << "evolution matrix: " << std::endl << U_t_t.slice((n_times - 1), 0) << std::endl;


    cx_cube occupancies(N, N, n_times, fill::zeros);
    for(int i = 0; i < n_times; ++i){
        occupancies.slice(i) = U_t_t.slice(i * coarsing_factor, 0) * fermi_function * U_t_t.slice(0, i * coarsing_factor);
    }

    mat local_occupancies(n_flavor, n_times, fill::zeros);
    for(int flavor = 0; flavor < n_flavor; ++flavor){
        cx_rowvec row = occupancies.tube(flavor * (n_k + 1), flavor * (n_k + 1));
        local_occupancies.row(flavor) = real(row);
    }

    local_occupancies.save(data.params_name + "_benchmark_occupancies.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_occupancies.txt", local_occupancies, {}, true);

    mat double_occupancy(1, n_times, fill::zeros);
    double_occupancy.row(0) = local_occupancies.row(0) % local_occupancies.row(1);
    double_occupancy.save(data.params_name + "_benchmark_double_occupancy.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_double_occupancy.txt", double_occupancy, {}, true);





    if (n_k % 2 == 0) {
        mat currents(n_flavor, n_times, fill::zeros);
        for (int flavor = 0; flavor < n_flavor; ++flavor) {
            for (int i = 0; i < n_times; ++i){
                for (int l = 0; l < n_k / 2; ++l){
                    currents(flavor, i) += -2 * std::imag(data.hybridizations[l][flavor][i]
                                           * occupancies(flavor * (n_k + 1) + 1 + l, flavor * (n_k + 1), i));
                }
            }
        }
        currents.save(data.params_name + "_benchmark_currents.out", raw_ascii);
        save_to_row::save(data.params_name + "_benchmark_currents.txt", currents, {}, true);
    }

    auto convert = [n_k, n_flavor](const cx_tetracube &input){
        int n_times_i = input.n_slices_1();
        int n_times_j = input.n_slices_2();
        cx_tetracube output{n_times_i, n_times_j, n_flavor, n_flavor};
        for(int a = 0; a < n_flavor; ++a){
            for(int b = 0; b < n_flavor; ++b){
                for(int i = 0; i < n_times_i; ++i){
                    for(int j = 0; j < n_times_j; ++j){
                        output(i, j, a, b) = input(a * (n_k + 1), b * (n_k + 1), i ,j);
                    }
                }
            }
        }
        return output;
    };

    cx_tetracube greater_gf(N, N, n_times, n_times);

    for(int i = 0; i < n_times; ++i){
        for(int j = 0; j < n_times; ++j){
            greater_gf.slice(i, j) = -I * U_t_t.slice(i * coarsing_factor, 0) * hole_fermi_function * U_t_t.slice(0, j * coarsing_factor);
        }
    }

    convert(greater_gf)().save(data.params_name + "_benchmark_greater_t_t_gf.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_gf_greater_t_t.txt", convert(greater_gf), {}, true);


    cx_tetracube lesser_gf(N, N, n_times, n_times);

    for(int i = 0; i < n_times; ++i){
        for(int j = 0; j < n_times; ++j){
            lesser_gf.slice(i, j) = I * U_t_t.slice(i * coarsing_factor, 0) * fermi_function * U_t_t.slice(0, j * coarsing_factor);
        }
    }

    convert(lesser_gf)().save(data.params_name + "_benchmark_lesser_t_t_gf.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_gf_lesser_t_t.txt", convert(lesser_gf), {}, true);

    std::cout << "Created real time GFs." << std::endl;

    cx_tetracube tau_t_gf(N, N, n_imag_times, n_times);
    cx_tetracube tau_gf(N, N, n_imag_times, 1);

    for(int i = 0; i < n_imag_times; ++i){
        for(int j = 0; j < n_times; ++j){
            tau_t_gf.slice(i, j) = -I * U_tau.slice(i * coarsing_factor) * hole_fermi_function * U_t_t.slice(0, j * coarsing_factor);
            if(j == 0){
                tau_gf.slice(i, 0) = -I * U_tau.slice(i * coarsing_factor) * hole_fermi_function;
            }
        }
    }

    convert(tau_t_gf)().save(data.params_name + "_benchmark_tau_t_gf.out", raw_ascii);
    convert(tau_gf)().save(data.params_name + "_benchmark_tau_gf.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_gf_tau_t.txt", convert(tau_t_gf), {}, true);
    save_to_row::save(data.params_name + "_benchmark_gf_tau.txt", convert(tau_gf), {}, true);



    std::cout << "Created mixing and imaginary time time GFs." << std::endl;

    const int n_matsubara = 100;

    cx_tetracube iw_t_gf(N, N, n_matsubara, n_times);

    for(int n = 0; n < n_matsubara; ++n){
        for(int j = 0; j < n_times; ++j){
            if(j == 0){
                iw_t_gf.slice(n, 0) = I * inv(I * (2 * n + 1.) * PI / data.beta * eye(size(h_0)) - h_0);
            }
            else{
                iw_t_gf.slice(n, j) = iw_t_gf.slice(n, 0) * U_t_t.slice(0, j * coarsing_factor);
            }
        }
    }

    convert(iw_t_gf)().save(data.params_name + "_benchmark_iw_t_gf.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_gf_iwn_t.txt", convert(iw_t_gf), {}, true);


    std::cout << "Created Matsubara GFs." << std::endl;

    const int n_w = 200;
    const double D = 5;
    const double dw = 2 * D / (n_w - 1);
    const double eta = 0.02;

    cx_tetracube spectral_function(N, N, n_w, 1);
    for(int n = 0; n < n_w; ++n){
        spectral_function.slice(n, 0) = conv_to<cx_mat>::from(
                -1. / PI * imag(inv((-D + n * dw + I * eta) * eye(size(h_0)) - h_0)));
    }

    convert(spectral_function)().save(data.params_name + "_benchmark_spectral_function.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_spectral_function.txt", convert(spectral_function), {}, true);

}

void InputGenerator::generate_interacting_benchmark_data() {

    std::cout << "Generating benchmark data (any U, n_bath_up = 1)" << std::endl;

    int n_flavor = data.local_energies.size();
    assert(n_flavor == 2);

    int n_k = data.bath_energies.size();
    if(n_k != 1){
        std::cout << "Input n_bath_up != 1, exiting" << std::endl;
        return;
    }

    std::cout << "Number of bath sites: " << n_k << std::endl;
    std::cout << "Number of flavors: " << n_flavor << std::endl;
    std::cout << "Number of time points: " << data.n_times << std::endl;

    // Generate h(t)
    int N = 16;
    cx_cube h(N, N, data.n_times, fill::zeros);
    const auto& eps_up = conv_to<cx_vec>::from(data.local_energies[0]);
    const auto& eps_do = conv_to<cx_vec>::from(data.local_energies[1]);
    const auto& eps_b_up = conv_to<cx_vec>::from(data.bath_energies[0][0]);
    const auto& eps_b_do = conv_to<cx_vec>::from(data.bath_energies[0][1]);
    const auto& V_up = conv_to<cx_vec>::from(data.hybridizations[0][0]);
    const auto& V_do = conv_to<cx_vec>::from(data.hybridizations[0][1]);
    const auto& U = conv_to<cx_vec>::from(data.Hubbard_U);

    for(int flavor = 0; flavor < n_flavor; ++flavor){
        h.tube(1, 1)= eps_up;
        h.tube(1, 3)= V_up;
        h.tube(2, 2)= eps_do;
        h.tube(2, 4)= V_do;
        h.tube(3, 3)= eps_b_up;
        h.tube(4, 4)= eps_b_do;
        h.tube(5, 5)= U + eps_do + eps_up;
        h.tube(5, 9)= V_do;
        h.tube(5, 10)= -V_up;
        h.tube(6, 6)= eps_b_do + eps_b_up;
        h.tube(6, 9)= V_up;
        h.tube(6, 10)= -V_do;
        h.tube(7, 7)= eps_b_up + eps_up;
        h.tube(8, 8)= eps_b_do + eps_do;
        h.tube(9, 9)= eps_b_do + eps_up;
        h.tube(10, 10)= eps_b_up + eps_do;
        h.tube(11, 11)= U + eps_b_up + eps_do + eps_up;
        h.tube(11, 13)= -V_do;
        h.tube(12, 12)= U + eps_b_do + eps_do + eps_up;
        h.tube(12, 14)= -V_up;
        h.tube(13, 13)= eps_b_do + eps_b_up + eps_up;
        h.tube(14, 14)= eps_b_do + eps_b_up + eps_do;
        h.tube(15, 15)= U + eps_b_do + eps_b_up + eps_do + eps_up;
    }

    for(int i = 0; i < h.n_slices; ++i){
        h.slice(i) = h.slice(i) + h.slice(i).t();
        h.slice(i).diag() /= 2.;
    }

    std::cout << "h_t: " << std::endl << h.slice(data.n_times - 1) << std::endl;


    EquationsOfMotion<> eom(h, data.dt, data.t_max, data.real_grid_size);
    eom.run();
    double error = eom.check_unitarity();
    std::cout << "Testing unitarity, maximal error: " << error << std::endl;

    auto U_t_t = eom.get_result();
    std::cout << U_t_t.slice(0, 0) << std::endl;

    cx_mat h_0(N, N, fill::zeros);
    h_0 = h.slice(0);

    cx_mat thermal_matrix = expmat_sym(-data.beta * h_0);
    double Z = trace(thermal_matrix).real();
    std::cout << "h_0: " << std::endl << h_0 << std::endl;
    std::cout << "fermi: " << std::endl << real(thermal_matrix / Z) << std::endl;

    int coarsing_factor = 1;
    int n_times = data.real_grid_size / coarsing_factor;

    std::cout << "evolution_matrix: " << std::endl << U_t_t.slice((n_times - 1), 0) << std::endl;


    cx_cube occupancies(N, N, n_times, fill::zeros);
    for(int i = 0; i < n_times; ++i){
        occupancies.slice(i) = U_t_t.slice(i * coarsing_factor, 0) * thermal_matrix / Z * U_t_t.slice(0, i * coarsing_factor);
    }

    mat local_occupancies(n_flavor, n_times, fill::zeros);
    mat double_occupancy(1, n_times, fill::zeros);
    mat currents(n_flavor, n_times, fill::zeros);
    cx_rowvec row;
    // Spin-up local level occ
    row.zeros(n_times);
    for(int state_index : {1, 5, 7, 9, 11, 12, 13, 15}){
        row += occupancies.tube(state_index, state_index);
    }
    local_occupancies.row(0) = real(row);
    // Spin-down local level occ
    row.zeros(n_times);
    for(int state_index : {2, 5, 8, 10, 11, 12, 14, 15}){
        row += occupancies.tube(state_index, state_index);
    }
    local_occupancies.row(1) = real(row);
    local_occupancies.save(data.params_name + "_benchmark_occupancies.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_occupancies.txt", local_occupancies, {}, true);


    // Double local level occ
    row.zeros(n_times);
    for(int state_index : {5, 11, 12, 15}){
        row += occupancies.tube(state_index, state_index);
    }
    double_occupancy.row(0) = real(row);
    double_occupancy.save(data.params_name + "_benchmark_double_occupancy.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_double_occupancy.txt", double_occupancy, {}, true);


    // Put hybridizations on the proper grid

    cx_vec hyb_up_old_grid = conv_to<cx_vec>::from(data.hybridizations[0][0]);
    cx_vec hyb_do_old_grid = conv_to<cx_vec>::from(data.hybridizations[0][1]);
    auto hyb_interp_up = Interpolation1D<cx_double>(hyb_up_old_grid, data.dt);
    auto hyb_interp_do = Interpolation1D<cx_double>(hyb_do_old_grid, data.dt);
    cx_rowvec hyb_up(n_times);
    cx_rowvec hyb_do(n_times);
    int dt = data.t_max / (n_times - 1);
    for(int i = 0; i < n_times; ++i){
        hyb_up(i) = hyb_interp_up(i * dt);
        hyb_do(i) = hyb_interp_do(i * dt);
    }

    // Current up
    row.zeros(n_times);
    row += occupancies.tube(3, 1);
    row += occupancies.tube(6, 9);
    row += -occupancies.tube(10, 5);
    row += -occupancies.tube(14, 12);

    currents.row(0) = -2 * imag(hyb_up % row);

    // Current down
    row.zeros(n_times);
    row += occupancies.tube(4, 2);
    row += -occupancies.tube(6, 10);
    row += occupancies.tube(9, 5);
    row += -occupancies.tube(13, 11);
    currents.row(1) = -2 * imag(hyb_do % row);

    currents.save(data.params_name + "_benchmark_currents.out", raw_ascii);
    save_to_row::save(data.params_name + "_benchmark_currents.txt", currents, {}, true);
}
