#include "BathWeight.h"

using namespace arma;

cx_double BathWeight::get_w_bath(){
    return get_sign_factor() * det_value;
}

cx_double BathWeight::get_det_ratio(){
    return det_ratio;
}

cx_double BathWeight::cut_hybridization_line(int i, int j) const {
    return M(j, i);
}

void BathWeight::generate_from_conf(){

    M.set_size(conf->get_k(), conf->get_k());
    int n_flavor = conf->get_n_flavor();

    int i = 0;
    int j = 0;
    for(int c_flavor = 0; c_flavor < n_flavor; ++c_flavor) {
        for(const auto& t_c : conf->get_c_times(c_flavor)) {
            j = 0;
            for(int a_flavor = 0; a_flavor < n_flavor; ++a_flavor) {
                for(const auto& t_a : conf->get_a_times(a_flavor)) {
                    M(i, j) = delta(c_flavor, a_flavor, t_c, t_a);
                    ++j;
                }
            }
            ++i;
        }
    }
    det_value = det(M);
//    if(std::abs(det_value) <= std::numeric_limits<double>::min())
//        throw std::runtime_error("Delta matrix is singular");
    M = inv(M);
}

void BathWeight::attempt_insertion(int flavor, const contour_time &t_c, const contour_time &t_a){
    col.set_size(conf->get_k());
    row.set_size(conf->get_k());
    int n_flavor = conf->get_n_flavor();

    int i = 0;
    for(int c_flavor = 0; c_flavor < n_flavor; ++c_flavor) {
        for(const auto &old_t_c : conf->get_c_times(c_flavor)) {
            col(i) = delta(c_flavor, flavor, old_t_c, t_a);
            i++;
        }
    }
    int j = 0;
    for(int a_flavor = 0; a_flavor < n_flavor; ++a_flavor) {
        for(const auto &old_t_a : conf->get_a_times(a_flavor)) {
            row(j) = delta(flavor, a_flavor, t_c, old_t_a);
            j++;
        }
    }

    scalar(0) = delta(flavor, flavor, t_c, t_a);
    scalar = scalar - row * M * col;
    // Note: this is correct only up to a external_sign! The external_sign is corrected in confirm_insertion method
    det_ratio = scalar(0);
}

void BathWeight::attempt_removal(int i, int j){
    // Note: this is correct only up to a external_sign! The external_sign is corrected in confirm_removal method
    det_ratio = M(j, i);
}

void BathWeight::attempt_shift_t_a(int a_flavor, int j, const contour_time &new_t_a){
    int n_flavor = conf->get_n_flavor();
    col.set_size(conf->get_k());

    int i = 0;
    for(int c_flavor = 0; c_flavor < n_flavor; ++c_flavor){
        for(const auto &t_c : conf->get_c_times(c_flavor)){
            col(i) = delta(c_flavor, a_flavor, t_c, new_t_a);
            i++;
        }
    }
    scalar = M.row(j) * col;
    // Note: this is correct only up to a external_sign! The external_sign is corrected in confirm_shift method
    det_ratio = scalar(0);
}


void BathWeight::attempt_shift_t_c(int c_flavor, int i, const contour_time &new_t_c){
    // NEW METHOD - NOT TESTED
    int n_flavor = conf->get_n_flavor();
    row.set_size(conf->get_k());

    int j = 0;
    for(int a_flavor = 0; a_flavor < n_flavor; ++a_flavor){
        for(const auto &t_a : conf->get_a_times(a_flavor)){
            row(j) = delta(c_flavor, a_flavor, new_t_c, t_a);
            j++;
        }
    }
    scalar = row * M.col(i);
    // Note: this is correct only up to a external_sign! The external_sign is corrected in confirm_shift method
    det_ratio = scalar(0);
}

void BathWeight::confirm_insertion(int i, int j){
    scalar = inv(scalar);
    col = -scalar(0) * M * col;
    row = -scalar(0) * row * M;
    M = M + col * row / scalar(0);
    M.insert_rows(j, row);
    col.insert_rows(j, scalar);
    M.insert_cols(i, col);
    // Note: external_sign correction here!
    det_ratio *= ((i + j) % 2 == 0 ? 1 : -1);
    det_value *= det_ratio;
}

void BathWeight::confirm_removal(int i, int j){
    col = M.col(i);
    M.shed_col(i);
    scalar(0) = col(j);
    col.shed_row(j);
    row = M.row(j);
    M.shed_row(j);
    M = M - col * row / scalar(0);
    // Note: external_sign correction here!
    det_ratio *= ((i + j) % 2 == 0 ? 1 : -1);
    det_value *= det_ratio;
}

void BathWeight::confirm_shift_t_a(int a_flavor, int j, int winding_direction){
    row = M.row(j) / scalar(0);
    M = M - M * col * row;
    M.row(j) = row;
    // If changing winding the first/last column must go to end/beginning since the a_times must be ordered
    int flavor_order = conf->get_flavor_order(a_flavor);
    if(winding_direction != 0 && flavor_order > 1){
        move_row(M, j, j - winding_direction * (flavor_order - 1));
        det_ratio *= ((flavor_order - 1) % 2 == 0 ? 1 : -1);
    }
    det_value *= det_ratio;
}

void BathWeight::confirm_shift_t_c(int c_flavor, int i, int winding_direction) {
    // NEW METHOD - NOT TESTED
    col = M.col(i) / scalar(0);
    M = M - col * row * M;
    M.col(i) = col;
    // If changing winding the first/last column must go to end/beginning since the a_times must be ordered
    int flavor_order = conf->get_flavor_order(c_flavor);
    if(winding_direction != 0 && flavor_order > 1){
        move_col(M, i, i - winding_direction * (flavor_order - 1));
        det_ratio *= ((flavor_order - 1) % 2 == 0 ? 1 : -1);
    }
    det_value *= det_ratio;
}

cx_double BathWeight::get_sign_factor() const{
    // external_sign factor = i^k
    switch(conf->get_k() % 4){
        case 0 : return 1.;
        case 1 : return I;
        case 2 : return -1.;
        case 3 : return -I;
    }
}

std::ostream &operator<<(std::ostream &os, const BathWeight &weight) {
    using std::endl;
    os << "*** Bath weight ***" << endl;
    os << "M: " << weight.M;
    os << "det_value: " << std::abs(weight.det_value) << " * e^i( " << std::arg(weight.det_value) / PI << " pi )" << endl;
    os << "sign_factor_phase: " << std::arg(weight.get_sign_factor()) / PI << " pi" << endl;
    return os;
}

