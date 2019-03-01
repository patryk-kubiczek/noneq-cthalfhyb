#include "DiscreteBathWeight.h"

void DiscreteBathWeight::set_n_blocks(int n_blocks) {
    discrete_bath = DiscreteBath(p->eps_0_up, p->eps_up, p->U_0, p->U,
                                 p->V_0_up, p->V_up, p->eps_bath_0_up, p->eps_bath_up,
                                 beta, t_max, n_blocks);
}

cx_long_double DiscreteBathWeight::calculate_new_trace() {
    //Compute the trace = phi * ... * phi * det
    int n_ini = conf->get_initial_n();
    auto times = conf->get_contour_times();
    new_trace_value = discrete_bath.calculate_z(times, n_ini);
    new_trace_value *= calculate_phi(times, n_ini);
    return new_trace_value;
}

cx_long_double DiscreteBathWeight::get_trace_ratio() const {
    return new_trace_value / trace_value;
}

void DiscreteBathWeight::finalize_update() {
    trace_value = new_trace_value;
    generate_permutation_sign();
}

cx_long_double DiscreteBathWeight::get_w_loc() const {
    return cx_long_double::value_type(permutation_sign) * trace_value;
}

arma::cx_mat DiscreteBathWeight::density_matrix_at_t_max(int *n_down) {
    return discrete_bath.calculate_density_matrix_at_t_max(conf->get_initial_n(), conf->get_contour_times(), n_down);
}

cx_double DiscreteBathWeight::calculate_phi(const std::vector<contour_time> &times, int n_ini) const {
    cx_double phi = 1.;
    contour_time ct_end;
    contour_time ct_start = {PLUS, 0};
    cx_double time_diff;
    for(const auto &ct : times){
        ct_end = ct;
        if(n_ini == 1) {
            time_diff = (ct_end.branch == IMAG ? -I : 1.) * ct_end.t - (ct_start.branch == IMAG ? -I : 1.) * ct_start.t;
            phi *= std::exp(-I * time_diff.real() * eps_down + time_diff.imag() * eps_0_down);
        }
        ct_start = ct;
        n_ini = (n_ini + 1) % 2;
    }
    if(n_ini == 1) {
        time_diff = -I * beta - (ct_start.branch == IMAG ? -I : 1.) * ct_start.t;
        phi *= std::exp(-I * time_diff.real() * eps_down + time_diff.imag() * eps_0_down);
    }
    return phi;
}

void DiscreteBathWeight::generate_from_conf() {
    trace_value = calculate_new_trace();
    generate_permutation_sign();
}

void DiscreteBathWeight::generate_permutation_sign() {
    if(conf->get_initial_n() == 1)
        permutation_sign = 1;
    else
        permutation_sign = (conf->get_k() % 2 == 0 ? 1 : -1);
}

std::ostream& DiscreteBathWeight::print(std::ostream &os) const {
    using std::endl;
    os << "*** Local Weight (Discrete Bath) ***" << endl;
    os << "trace_value: " << std::abs(trace_value) << " * e^i( " << std::arg(trace_value) / PI << " pi )" << endl;
    os << "permutation_sign: " << permutation_sign << endl;
    return os;
}


