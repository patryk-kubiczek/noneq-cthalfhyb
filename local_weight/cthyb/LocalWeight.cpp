#include "LocalWeight.h"

using namespace arma;

LocalWeight::LocalWeight(const Model &m, const Configuration &conf, local_weight_type type)
        : LocalWeightBase(type),
          t_max(m.t_max), beta(m.beta),
          p(ContourPropagatorWrapper::create(m, type)),
          c_operators(&m.c_operator_matrices), a_operators(&m.a_operator_matrices),
          conf(&conf) {
    auto n_flavor = conf.get_n_flavor();
    c_matrix_iters.resize(n_flavor);
    a_matrix_iters.resize(n_flavor);
    generate_from_conf();
}

void LocalWeight::generate_from_conf() {
    // Import operators from conf and sort them increasingly w.r.t. contour time
    struct ContourOperator {
        operator_type type;
        int flavor;
        int global_index;
        int relative_index;
        contour_time time;
    };
    std::vector<ContourOperator> contour_operators;
    auto n_flavor = conf->get_n_flavor();
    int index = 0;
    for (int flavor = 0; flavor < n_flavor; ++flavor) {
        // Number of creation and anihiliation operators of the same flavor must be equal
        for (int i = 0; i < conf->get_flavor_order(flavor); ++i) {
            contour_operators.push_back(
                    {annihilation, flavor, index++, i, conf->get_contour_time(annihilation, flavor, i)});
            contour_operators.push_back({creation, flavor, index++, i, conf->get_contour_time(creation, flavor, i)});
        }
    }
    std::sort(begin(contour_operators), end(contour_operators),
              [](const auto &lhs, const auto &rhs) { return lhs.time < rhs.time; });

    // Calculate the permutation external_sign
    permutation_sign = get_permutation_sign([&contour_operators](int n) { return contour_operators[n].global_index; },
                                            contour_operators.size());

    // Generate local_matrices and store pointers to fermionic operators
    local_matrices.clear();
    c_matrix_iters.clear();
    a_matrix_iters.clear();
    for (int flavor = 0; flavor < n_flavor; ++flavor) {
        c_matrix_iters.emplace_back(matrix_iters_t(conf->get_flavor_order(flavor)));
        a_matrix_iters.emplace_back(matrix_iters_t(conf->get_flavor_order(flavor)));
    }

    contour_time start_time = {PLUS, 0};
    for (auto &op : contour_operators) {
        local_matrices.emplace_back(make_local_propagator(op.time, start_time));
        local_matrices.emplace_back(make_local_operator(op.type, op.flavor, op.time));
        (op.type == annihilation ? a_matrix_iters : c_matrix_iters)[op.flavor][op.relative_index] = prev(
                end(local_matrices));
        start_time = op.time;
    }
    local_matrices.emplace_back(make_local_propagator({IMAG, beta}, start_time));

    trace_value = calculate_new_trace();
}

cx_long_double LocalWeight::calculate_new_trace() {
    int matrix_element = (is_zero_order_line_present() ? get_contributing_matrix_element() : -1);
    new_trace_value = trace_of_propagators(matrix_element);
    return new_trace_value;
}

cx_long_double LocalWeight::get_trace_ratio() const {
    return new_trace_value / trace_value;
}

cx_long_double LocalWeight::get_w_loc() const {
    return cx_long_double::value_type(permutation_sign) * trace_value;
}

void LocalWeight::attempt_add_operator(operator_type type, int flavor, const contour_time &time) {
    // This function updates local_matrices after an operator is added

    // Find the iterator to the propagator we wish to remove
    auto mat_iter = std::upper_bound(begin(local_matrices), end(local_matrices), time,
                                     [](const auto &t, const auto &mat) { return t < mat->get_end_time(); });
    const auto start_time = (*mat_iter)->get_start_time();
    const auto end_time = (*mat_iter)->get_end_time();

    auto removed_p_iter = mat_iter;
    ++mat_iter;

    // Remove the propagator by moving it to another list
    export_to_vector_of_lists(old_local_matrices, local_matrices, removed_p_iter);

    // Insert a new propagator, the new operator and another new propagator
    mat_iter = local_matrices.insert(mat_iter, make_local_propagator(time, start_time));
    auto first_new_mat_iter = mat_iter;
    ++mat_iter;
    mat_iter = local_matrices.insert(mat_iter, make_local_operator(type, flavor, time));
    ++mat_iter;
    mat_iter = local_matrices.insert(mat_iter, make_local_propagator(end_time, time));
    new_mat_iters.emplace_back(make_pair(first_new_mat_iter, next(mat_iter)));

    // Returns iterator to the new operator
    (type == annihilation ? temporary_a_iter : temporary_c_iter) = next(first_new_mat_iter);
}

void LocalWeight::attempt_remove_operator(operator_type op_type, int flavor, int relative_index) {
    // This function updates local_matrices after an operator is removed

    // Find the iterator to the operator we wish to remove
    auto removed_op_iter = (op_type == annihilation ? a_matrix_iters : c_matrix_iters)[flavor][relative_index];
    const auto start_time = (*prev(removed_op_iter))->get_start_time();
    const auto end_time = (*next(removed_op_iter))->get_end_time();
    auto mat_iter = next(next(removed_op_iter));

    // Remove the first propagator, the operator and the second propagator by moving them to another list
    export_to_vector_of_lists(old_local_matrices, local_matrices, prev(removed_op_iter), mat_iter);

    // Insert a new propagator
    mat_iter = local_matrices.insert(mat_iter, make_local_propagator(end_time, start_time));

    new_mat_iters.emplace_back(make_pair(mat_iter, next(mat_iter)));
}

void LocalWeight::confirm_add_operator(operator_type op_type, int flavor, int relative_index) {
    const auto &op_iter = (op_type == annihilation ? temporary_a_iter : temporary_c_iter);
    // Should be executed after updating conf
    insert_to_vector((op_type == annihilation ? a_matrix_iters : c_matrix_iters)[flavor], relative_index, op_iter);
}

void LocalWeight::confirm_remove_operator(operator_type op_type, int flavor, int relative_index) {
    // Should be executed after updating conf
    erase_from_vector((op_type == annihilation ? a_matrix_iters : c_matrix_iters)[flavor], relative_index);
}

void LocalWeight::finalize_update() {
    trace_value = new_trace_value;
    generate_permutation_sign();
    old_local_matrices.clear();
    new_mat_iters.clear();
}

void LocalWeight::cancel_operations() {
    local_matrices_t::iterator mat_iter;

    assert(old_local_matrices.size() == new_mat_iters.size());
    for (int r = old_local_matrices.size() - 1; r >= 0; --r) {
        mat_iter = local_matrices.erase(new_mat_iters[r].first, new_mat_iters[r].second);
        local_matrices.splice(mat_iter, old_local_matrices[r]);
    }
    old_local_matrices.clear();
    new_mat_iters.clear();
}

std::unique_ptr<LocalWeight::LocalFermionicOperator>
LocalWeight::make_local_operator(operator_type type, int flavor, const contour_time &time) {
    return std::make_unique<LocalFermionicOperator>(time,
                                                    (type == annihilation ? *a_operators : *c_operators).slice(flavor));
}

std::unique_ptr<LocalWeight::LocalPropagator>
LocalWeight::make_local_propagator(const contour_time &end_time, const contour_time &start_time) {
    return std::make_unique<LocalPropagator>(end_time, start_time, p(end_time, start_time));
}

bool LocalWeight::is_zero_order_line_present() const {
    for (int flavor = 0; flavor < conf->get_n_flavor(); ++flavor) {
        if (conf->get_flavor_order(flavor) == 0) return true;
    }
    return false;
}

int LocalWeight::get_contributing_matrix_element() const {
    int matrix_element = 0;
    int power_of_two = 1;
    for (int flavor = 0; flavor < conf->get_n_flavor(); ++flavor) {
        matrix_element += power_of_two * conf->get_initial_n(flavor);
        power_of_two *= 2;
    }
    return matrix_element;
}

void LocalWeight::generate_permutation_sign() {
    // Import operators from conf and sort them increasingly w.r.t. contour time
    struct ContourOperator {
        int global_index;
        contour_time time;
    };
    std::vector<ContourOperator> contour_operators;
    auto n_flavor = conf->get_n_flavor();
    int index = 0;
    for (int flavor = 0; flavor < n_flavor; ++flavor) {
        // Number of creation and anihiliation operators of the same flavor must be equal
        for (int i = 0; i < conf->get_flavor_order(flavor); ++i) {
            assert(conf->get_global_index(annihilation, flavor, i) == index);
            contour_operators.push_back({index++, conf->get_contour_time(annihilation, flavor, i)});
            assert(conf->get_global_index(creation, flavor, i) == index);
            contour_operators.push_back({index++, conf->get_contour_time(creation, flavor, i)});
        }
    }
    std::sort(begin(contour_operators), end(contour_operators),
              [](const auto &lhs, const auto &rhs) { return lhs.time < rhs.time; });

    // Calculate the permutation external_sign
    permutation_sign = get_permutation_sign([&contour_operators](int n) { return contour_operators[n].global_index; },
                                            contour_operators.size());
}

cx_double LocalWeight::trace_of_propagators(int matrix_element = -1) {
    const auto N = local_matrices.front()->get().n_rows;
    cx_mat product(N, N, fill::eye);
    for (const auto &m : local_matrices) {
        product = m->get() * product;
    }
    if (matrix_element == -1) return trace(diagmat(product));
    else return product(matrix_element, matrix_element);
}

std::ostream &LocalWeight::print(std::ostream &os) const {
    using std::endl;
    os << "*** Local Weight ***" << endl;
    os << "local matrices: [";
    bool first = true;
    for (const auto &m : local_matrices) {
        if (!first)
            os << ", ";
        else
            first = false;
        os << "(" << m->get_start_time().s(t_max) << ", " << m->get_end_time().s(t_max) << ")";
    }
    os << "]" << endl;
    os << "number of local matrices: " << local_matrices.size() << endl;
    os << "number of a matrix iters: " << endl;
    int n_flavor = conf->get_n_flavor();
    for (int a = 0; a < n_flavor; ++a) {
        os << " flavor no. " << a << ": " << a_matrix_iters[a].size() << endl;
    }
    os << "number of c matrix iters: " << endl;
    for (int a = 0; a < n_flavor; ++a) {
        os << " flavor no. " << a << ": " << c_matrix_iters[a].size() << endl;
    }
    os << "trace_value: " << std::abs(trace_value) << " * e^i( "
       << std::arg(trace_value) / PI << " pi )" << endl;
    os << "permutation_sign: " << permutation_sign << endl;
    os << "number_of_old_matrices: " << old_local_matrices.size() << endl;
    os << "number_of_new_mat_iters: " << new_mat_iters.size();
    return os;
}



