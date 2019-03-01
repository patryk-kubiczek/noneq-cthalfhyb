#include "LocalWeight.h"

using namespace arma;

LocalWeight::LocalWeight(const Model &m, const Configuration &conf, local_weight_type type)
        : LocalWeightBase(type),
          t_max(m.t_max), beta(m.beta),
          u_0(ContourUWrapper<0>::create(m, type)), u_1(ContourUWrapper<1>::create(m, type)),
          phi_1(ContourPhi1Wrapper::create(m, type)),
          conf(&conf) {
    generate_from_conf();
}

void LocalWeight::generate_from_conf() {
    assert(conf->get_n_flavor() == 1);

    // Import operators from conf and sort them increasingly w.r.t. contour time
    struct ContourOperator {
        operator_type type;
        int relative_index;
        contour_time time;
    };
    std::vector<ContourOperator> contour_operators;

    for (int i = 0; i < conf->get_k(); ++i) {
        contour_operators.push_back({annihilation, i, conf->get_contour_time(annihilation, 0, i)});
        contour_operators.push_back({creation, i, conf->get_contour_time(creation, 0, i)});
    }

    std::sort(begin(contour_operators), end(contour_operators),
              [](const auto &lhs, const auto &rhs) { return lhs.time < rhs.time; });

    // Calculate the permutation sign
    generate_permutation_sign();

    // Generate evolution_matrices and store pointers to fermionic operators
    evolution_matrices.clear();
    c_matrix_iters.resize(conf->get_k());
    a_matrix_iters.resize(conf->get_k());

    contour_time start_time = {PLUS, 0};
    for (auto &op : contour_operators) {
        if (op.type == annihilation) {
            evolution_matrices.emplace_back(make_evolution_matrix(1, op.time, start_time));
            a_matrix_iters[op.relative_index] = prev(end(evolution_matrices));
        }
        else {
            evolution_matrices.emplace_back(make_evolution_matrix(0, op.time, start_time));
            c_matrix_iters[op.relative_index] = prev(end(evolution_matrices));
        }
        start_time = op.time;
    }
    evolution_matrices.emplace_back(make_evolution_matrix(conf->get_initial_n(), {IMAG, beta}, start_time));

    //Compute the trace
    trace_value = calculate_new_trace();
}

cx_long_double LocalWeight::calculate_new_trace() {
    //Compute the trace = phi * ... * phi * det
    new_trace_value = determinant();
    for(const auto &mat : evolution_matrices)
        new_trace_value *= mat.get_phi();
    return new_trace_value;
}

cx_long_double LocalWeight::get_trace_ratio() const {
    return new_trace_value / trace_value;
}

cx_long_double LocalWeight::get_w_loc() const {
    return cx_long_double::value_type(permutation_sign) * trace_value;
}

cx_mat LocalWeight::density_matrix_at_t_max(int *n_down_at_t_max) {
    int n_down;
    std::vector<EvolutionMatrix> reordered_matrices;
    const int k = conf->get_k();
    if(k > 0) {
        auto segments = conf->get_segments();
        int i = 0;
        while(!segments[i].contains({PLUS, t_max})) i = (i + 1) % (2 * k);
        n_down = segments[i].n();
        reordered_matrices.emplace_back(make_evolution_matrix(n_down, segments[i].end_time(), {MINUS, t_max}));
        for(int j = (i + 1) % (2 * k); j != i; j = (j + 1) % (2 * k)){
            reordered_matrices.emplace_back(
                    make_evolution_matrix(segments[j].n(), segments[j].end_time(), segments[j].start_time()));
        }
        reordered_matrices.emplace_back(make_evolution_matrix(n_down, {PLUS, t_max}, segments[i].start_time()));
    }
    else {
        n_down = conf->get_initial_n();
        reordered_matrices.emplace_back(make_evolution_matrix(n_down, {IMAG, beta}, {MINUS, t_max}));
        reordered_matrices.emplace_back(make_evolution_matrix(n_down, {PLUS, t_max}, {PLUS, 0}));
    }

    const cx_mat one(size(reordered_matrices[0].get()), fill::eye);
    cx_mat U(one);
    for(const auto& matrix : reordered_matrices){
        U = matrix.get() * U;
    }
    if(n_down_at_t_max != nullptr)
        *n_down_at_t_max = n_down;
    return inv(one + U);
}


void LocalWeight::attempt_add_pair(int n_down, const contour_time &end_time, const contour_time &start_time) {
    // This function updates evolution_matrices after an operator pair is added

    // This function DOES NOT work for adding wound operator_pairs
    assert(start_time < end_time);

    // Find the iterator to the evolution matrix we wish to split
    auto mat_iter = std::upper_bound(begin(evolution_matrices), end(evolution_matrices), end_time,
                                     [](const auto &t, const auto &mat) { return t < mat.get_end_time(); });

    const auto old_start_time = (*mat_iter).get_start_time();
    const auto old_end_time = (*mat_iter).get_end_time();
    assert(start_time > old_start_time);

    auto removed_mat_iter = mat_iter;
    ++mat_iter;

    // Remove the evolution matrix by moving it to another list

    export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, removed_mat_iter);

    // Insert three new evolution matrices

    int old_n_down = (n_down + 1) % 2;
    mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(old_n_down, start_time, old_start_time));
    auto first_new_mat_iter = mat_iter;
    ++mat_iter;
    mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(n_down, end_time, start_time));
    ++mat_iter;
    mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(old_n_down, old_end_time, end_time));

    new_mat_iters.emplace_back(make_pair(first_new_mat_iter, next(mat_iter)));

    // Returns iterator to the middle evolution matrix
    temporary_iter_addition = next(first_new_mat_iter);
}

void LocalWeight::attempt_remove_pair(int n_down, int relative_index) {
    // This function updates evolution_matrices after an operator pair is removed

    // Find the iterator to the segment/antisegment evolution matrix we wish to remove
    auto removed_mat_iter = (n_down == 0 ? c_matrix_iters : a_matrix_iters)[relative_index];

    // This function DOES NOT work for removing wound operator pairs
    assert(removed_mat_iter->get_start_time() < removed_mat_iter->get_end_time());

    const auto start_time = (*prev(removed_mat_iter)).get_start_time();
    const auto end_time = (*next(removed_mat_iter)).get_end_time();
    auto mat_iter = next(next(removed_mat_iter));

    // Remove the previous evolution matrix, the evolution matrix corresponding to removed pair and the next evolution matrix by moving them to another list

    export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, prev(removed_mat_iter), mat_iter);

    // Insert a new evolution matrix
    mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix((n_down + 1) % 2, end_time, start_time));

    new_mat_iters.emplace_back(make_pair(mat_iter, next(mat_iter)));

    // Returns iterator to the new evolution matrix
    temporary_iter_removal = mat_iter;
}

void LocalWeight::attempt_add_wound_pair(int n_down, const contour_time &end_time, const contour_time &start_time) {
    // This function updates evolution_matrices after a wound operator pair is added
    assert(start_time > end_time);

    auto mat_iter = begin(evolution_matrices);

    if (get_current_order() > 0) {
        // Remove the first evolution matrix by moving it to another list
        auto old_time = (*mat_iter).get_end_time();
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter);

        // Insert two new evolution matrices at the beginning
        mat_iter = begin(evolution_matrices);
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(n_down, end_time, {PLUS, 0}));
        ++mat_iter;
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix((n_down + 1) % 2, old_time, end_time));

        // Insert the iterators to the new evolution matrices
        new_mat_iters.emplace_back(make_pair(prev(mat_iter), (next(mat_iter))));

        // Remove the last evolution matrix by moving it to another list
        mat_iter = prev(end(evolution_matrices));
        old_time = (*mat_iter).get_start_time();
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter);

        // Insert two new evolution matrices at the end
        mat_iter = end(evolution_matrices);
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix((n_down + 1) % 2, start_time, old_time));
        ++mat_iter;
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(n_down, {IMAG, beta}, start_time));

        // Insert the iterators to the new evolution matrices
        new_mat_iters.emplace_back(make_pair(prev(mat_iter), next(mat_iter)));
    }
    else {
        // Case in which a wound pair is added to an empty line needs separate treatment

        // Remove the full evolution matrix by moving it to another list
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter);

        // Insert three new evolution matrices
        mat_iter = begin(evolution_matrices);
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(n_down, end_time, {PLUS, 0}));
        ++mat_iter;
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix((n_down + 1) % 2, start_time, end_time));
        ++mat_iter;
        mat_iter = evolution_matrices.insert(mat_iter, make_evolution_matrix(n_down, {IMAG, beta}, start_time));

        // Insert the iterators to the new evolution matrices
        new_mat_iters.emplace_back(make_pair(prev(prev(mat_iter)), next(mat_iter)));
    }
}

void LocalWeight::attempt_remove_wound_pair(int n_down) {
    // This function updates evolution_matrices after a wound operator pair is removed
    assert(get_current_order() > 0);

    // Find the iterator to the first evolution matrix we wish to remove
    auto mat_iter = begin(evolution_matrices);

    if (get_current_order() > 1) {
        auto old_time = (*next(mat_iter)).get_end_time();
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter, next(next(mat_iter)));
        mat_iter = evolution_matrices.insert(begin(evolution_matrices),
                                             make_evolution_matrix((n_down + 1) % 2, old_time, {PLUS, 0}));
        new_mat_iters.emplace_back(make_pair(mat_iter, next(mat_iter)));

        mat_iter = prev(prev(end(evolution_matrices)));
        old_time = (*mat_iter).get_start_time();
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter, next(next(mat_iter)));
        mat_iter = evolution_matrices.insert(end(evolution_matrices),
                                             make_evolution_matrix((n_down + 1) % 2, {IMAG, beta}, old_time));
        new_mat_iters.emplace_back(make_pair(mat_iter, next(mat_iter)));
    }
    else {
        export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, mat_iter, end(evolution_matrices));
        mat_iter = evolution_matrices.insert(begin(evolution_matrices),
                                             make_evolution_matrix((n_down + 1) % 2, {IMAG, beta}, {PLUS, 0}));
        new_mat_iters.emplace_back(make_pair(mat_iter, next(mat_iter)));
    }
}

void LocalWeight::attempt_switch_n(int new_n_down) {
    assert(get_current_order() == 0);
    export_to_vector_of_lists(old_evolution_matrices, evolution_matrices, begin(evolution_matrices), end(evolution_matrices));
    evolution_matrices.insert(begin(evolution_matrices), make_evolution_matrix(new_n_down, {IMAG, beta}, {PLUS, 0}));
    new_mat_iters.emplace_back(make_pair(begin(evolution_matrices), end(evolution_matrices)));
}

void LocalWeight::confirm_add_pair(int n_down, int c_relative_index, int a_relative_index) {
    const auto& mat_iter = temporary_iter_addition;
    if (n_down == 0) {
        insert_to_vector(c_matrix_iters, c_relative_index, mat_iter);
        if (a_relative_index < a_matrix_iters.size())
            a_matrix_iters[a_relative_index] = next(mat_iter);
        insert_to_vector(a_matrix_iters, a_relative_index, prev(mat_iter));
    } 
    else {
        insert_to_vector(a_matrix_iters, a_relative_index, mat_iter);
        if (c_relative_index < c_matrix_iters.size())
            c_matrix_iters[c_relative_index] = next(mat_iter);
        insert_to_vector(c_matrix_iters, c_relative_index, prev(mat_iter));
    }
}

void LocalWeight::confirm_remove_pair(int n_down, int c_relative_index, int a_relative_index) {
    const auto& mat_iter = temporary_iter_removal;
    erase_from_vector(c_matrix_iters, c_relative_index);
    erase_from_vector(a_matrix_iters, a_relative_index);

    if (n_down == 0) {
        if (a_relative_index < a_matrix_iters.size())
            a_matrix_iters[a_relative_index] = mat_iter;
    }
    else {
        if (c_relative_index < c_matrix_iters.size())
            c_matrix_iters[c_relative_index] = mat_iter;
    }
}

void LocalWeight::confirm_add_wound_pair(int n_down) {
    if (n_down == 0) {
        insert_to_vector(c_matrix_iters, 0, begin(evolution_matrices));
        if (a_matrix_iters.size() > 0)
            a_matrix_iters[0] = next(begin(evolution_matrices));
        insert_to_vector(a_matrix_iters, a_matrix_iters.size(), prev(prev(end(evolution_matrices))));
    }
    else {
        insert_to_vector(a_matrix_iters, 0, begin(evolution_matrices));
        if (c_matrix_iters.size() > 0)
            c_matrix_iters[0] = next(begin(evolution_matrices));
        insert_to_vector(c_matrix_iters, c_matrix_iters.size(), prev(prev(end(evolution_matrices))));
    }
}

void LocalWeight::confirm_remove_wound_pair(int n_down) {
    if (n_down == 0) {
        erase_from_vector(c_matrix_iters, 0);
        erase_from_vector(a_matrix_iters, a_matrix_iters.size() - 1);
        if (a_matrix_iters.size() > 0)
            a_matrix_iters[0] = begin(evolution_matrices);
    }
    else {
        erase_from_vector(a_matrix_iters, 0);
        erase_from_vector(c_matrix_iters, c_matrix_iters.size() - 1);
        if (c_matrix_iters.size() > 0)
            c_matrix_iters[0] = begin(evolution_matrices);
    }
}


void LocalWeight::finalize_update() {
    generate_permutation_sign();
    trace_value = new_trace_value;
    old_evolution_matrices.clear();
    new_mat_iters.clear();
}

void LocalWeight::cancel_operations() {
    evolution_matrices_t::iterator mat_iter;

    assert(old_evolution_matrices.size() == new_mat_iters.size());
    for (int r = old_evolution_matrices.size() - 1; r >= 0; --r) {
        mat_iter = evolution_matrices.erase(new_mat_iters[r].first, new_mat_iters[r].second);
        evolution_matrices.splice(mat_iter, old_evolution_matrices[r]);
    }
    old_evolution_matrices.clear();
    new_mat_iters.clear();
}

LocalWeight::EvolutionMatrix
LocalWeight::make_evolution_matrix(int n_down, const contour_time &end_time, const contour_time &start_time) const {
    if(n_down == 0)
        return EvolutionMatrix{end_time, start_time, u_0(end_time, start_time), 1.};
    else
        return EvolutionMatrix{end_time, start_time, u_1(end_time, start_time), phi_1(end_time, start_time)};
}

int LocalWeight::get_current_order() const {
    return (evolution_matrices.size() - 1) / 2;
}

void LocalWeight::generate_permutation_sign() {
    if (conf->get_initial_n() == 1)
        permutation_sign = 1;
    else
        permutation_sign = (conf->get_k() % 2 == 0 ? 1 : -1);
}

cx_long_double LocalWeight::determinant() const {
    const mat one(size(evolution_matrices.front().get()), fill::eye);
    cx_mat product(size(one), fill::eye);
    for(const auto &m : evolution_matrices)
        product = m.get() * product;
    return std::exp(cx_long_double{log_det(one + product)});
}

std::ostream &LocalWeight::print(std::ostream &os) const {
    using std::endl;
    os << "*** Local Weight ***" << endl;
    os << "evolution matrices: [";
    bool first = true;
    for(const auto &m : evolution_matrices){
        if(!first)
            os << ", ";
        else
            first = false;
        os << "(" << m.get_start_time().s(t_max) << ", " << m.get_end_time().s(t_max) << ")" ;
    }
    os << "]" << endl;
    os << "number of evolution matrices: " << evolution_matrices.size() << endl;
    os << "number of old evolution matrices: " << old_evolution_matrices.size() << endl;
    os << "number of new matrices iterators: " << new_mat_iters.size() << endl;
    os << "number of c matrix iters: " << c_matrix_iters.size() << endl;
    os << "number of a matrix iters: " << a_matrix_iters.size() << endl;
    os << "trace_value: " << std::abs(trace_value) << " * e^i( "
       << std::arg(trace_value) / PI << " pi )" << endl;
    os << "permutation_sign: " << permutation_sign;
    return os;
}


