#include "Configuration.h"

Configuration::Configuration(double t_max, double beta, int n_flavor) : t_max(t_max), beta(beta) {
    c_times.resize(n_flavor);
    a_times.resize(n_flavor);
    initial_n.resize(n_flavor);
    flavor_order.resize(n_flavor);
}

//Configuration::Configuration(const Configuration &conf) : t_max(conf.t_max), beta(conf.beta) {
//    // Copy constructor only sets constant parts of conf
//    c_times.resize(conf.get_n_flavor());
//    a_times.resize(conf.get_n_flavor());
//    initial_n.resize(conf.get_n_flavor());
//    flavor_order.resize(conf.get_n_flavor());
//}
//
//Configuration &Configuration::operator=(Configuration &&conf) {
//    // Assignment sets the varying parts of conf
//    k = conf.k;
//    abs_weight = conf.abs_weight;
//    external_sign = conf.external_sign;
//    internal_sign = conf.internal_sign;
//    c_times = std::move(conf.c_times);
//    a_times = std::move(conf.a_times);
//    initial_n = std::move(conf.initial_n);
//    flavor_order = std::move(conf.flavor_order);
//    return *this;
//}
//
//Configuration &Configuration::operator=(const Configuration &conf) {
//    // Assignment sets the varying parts of conf
//    k = conf.k;
//    abs_weight = conf.abs_weight;
//    external_sign = conf.external_sign;
//    internal_sign = conf.internal_sign;
//    c_times = conf.c_times;
//    a_times = conf.a_times;
//    initial_n = conf.initial_n;
//    flavor_order = conf.flavor_order;
//    return *this;
//}

int Configuration::get_absolute_index(operator_type op_type, int flavor, int relative_index) const {
    int index = 0;
    for(int a = 0; a < flavor; ++a){
        assert(flavor_order[a] == c_times[a].size() && flavor_order[a] == a_times[a].size());
        index += flavor_order[a];
    }
    index += relative_index; // now we have the correct index among only operators of op_type
    return index;
}

int Configuration::get_global_index(operator_type op_type, int flavor, int relative_index) const {
    int index = get_absolute_index(op_type, flavor, relative_index);
    index *= 2; // now we include the other type of operators by doubling the index
    index += (op_type == annihilation ? 0 : 1); // creation operator is always one position further than anihilation
    return index;
}

int Configuration::get_relative_index(int global_index, operator_type op_type, int &flavor) const {
    op_type = (global_index % 2 == 0 ? annihilation : creation);
    global_index -= (op_type == annihilation ? 0 : 1);
    global_index /= 2;
    int a = 0;
    int new_index = global_index;
    while(true){
        global_index = new_index;
        new_index -= (op_type == annihilation ? a_times : c_times)[a].size();
        if(new_index < 0)
            break;
        ++a;
    }
    flavor = a;
    return global_index;
}

void Configuration::calculate_internal_sign() {
    //(-1)^k
    internal_sign = (k % 2 == 0 ? 1 : -1);
    //\lambda(t)
    int n_flavor = get_n_flavor();
    for(int flavor = 0; flavor < n_flavor; ++flavor){
        for(const auto& t : c_times[flavor])
            internal_sign *= t.integration_sign();
        for(const auto& t : a_times[flavor])
            internal_sign *= t.integration_sign();
    }
}



void Configuration::randomly_generate_general(int k, RandomNumberGenerator &rng) {
    this->k = k;
    int n_flavor = get_n_flavor();
    int flavor;
    for(flavor = 0; flavor < n_flavor; ++flavor){
        c_times[flavor].clear();
        a_times[flavor].clear();
    }
    for(int i = 0; i < k; ++i) {
        flavor = rng.random_int(0, n_flavor - 1);
        ++flavor_order[flavor];
        c_times[flavor].emplace_back(rng.random_contour_time({PLUS, 0}, {IMAG, beta}, t_max, beta));
        a_times[flavor].emplace_back(rng.random_contour_time({PLUS, 0}, {IMAG, beta}, t_max, beta));
    }
    for(flavor = 0; flavor < n_flavor; ++flavor){
        std::sort(begin(c_times[flavor]), end(c_times[flavor]));
        std::sort(begin(a_times[flavor]), end(a_times[flavor]));
        if(get_flavor_order(flavor) == 0) initial_n[flavor] = rng.random_int(0, 1);
        else initial_n[flavor] = get_initial_n(flavor);
    }
    calculate_internal_sign();
}

void Configuration::randomly_generate_segment_like(int k, RandomNumberGenerator &rng) {
    this->k = k;
    int n_flavor = get_n_flavor();
    int flavor;
    for(flavor = 0; flavor < n_flavor; ++flavor){
        c_times[flavor].clear();
        a_times[flavor].clear();
    }
    // First generate random times for each flavor not assigning them yet to any operators
    std::vector<std::vector<contour_time>> random_times(n_flavor);
    for(int i = 0; i < k; ++i) {
        flavor = rng.random_int(0, n_flavor - 1);
        ++flavor_order[flavor];
        random_times[flavor].emplace_back(rng.random_contour_time({PLUS, 0}, {IMAG, beta}, t_max, beta));
        random_times[flavor].emplace_back(rng.random_contour_time({PLUS, 0}, {IMAG, beta}, t_max, beta));
    }
    // Sort the generated times and interchangeably insert creation and annihilation times
    int n_on_contour;
    for(flavor = 0; flavor < n_flavor; ++flavor){
        std::sort(begin(random_times[flavor]), end(random_times[flavor]));
        // Randomly generate initial n on line
        n_on_contour = rng.random_int(0, 1);
        for(const auto& time : random_times[flavor]){
            if(n_on_contour == 0)
                c_times[flavor].push_back(time);
            else
                a_times[flavor].push_back(time);
            n_on_contour = (n_on_contour + 1) % 2;
        }
        if(get_flavor_order(flavor) == 0) initial_n[flavor] = rng.random_int(0, 1);
        else initial_n[flavor] = get_initial_n(flavor);
    }
    calculate_internal_sign();

}

void Configuration::randomly_generate(int k, RandomNumberGenerator &rng, bool segment_like) {
    // For OneSpin segment_like HAS TO be true
    randomly_generate_segment_like(k, rng);
}

std::vector<Configuration::Segment> Configuration::get_segments(int flavor) const {
    const int k = flavor_order[flavor];
    if(k == 0) return std::vector<Segment>{};
    const int n = (c_times[flavor][0] < a_times[flavor][0] ? 0 : 1);
    std::vector<Segment> segments;
    for(int i = 0; i < k; ++i){
        // argument order = (c_position, a_position, t_c, t_a, n, flavor)
        segments.emplace_back(Segment(i, i, c_times[flavor][i], a_times[flavor][i], (n + 1) % 2, flavor));
        segments.emplace_back(Segment((i + 1 - n) % k, (i + n) % k,
                                      c_times[flavor][(i + 1 - n) % k], a_times[flavor][(i + n) % k], n, flavor));
    }
    return segments;
}

std::vector<contour_time> Configuration::get_contour_times(int flavor) const {
    const int k = flavor_order[flavor];
    if(k == 0) return std::vector<contour_time>{};
    std::vector<contour_time> contour_times;
    const int n_ini = get_initial_n(flavor);
    for(int i = 0; i < k; ++i){
        contour_times.push_back((n_ini == 0 ? c_times : a_times)[flavor][i]);
        contour_times.push_back((n_ini == 0 ? a_times : c_times)[flavor][i]);
    }
    return contour_times;
}

int Configuration::get_initial_n(int flavor) const{
    if(flavor_order[flavor] > 0) return c_times[flavor][0] < a_times[flavor][0] ? 0 : 1;
    else return initial_n[flavor];
}


bool Configuration::has_operators_on_real_branch() const {
    for(int flavor = 0; flavor < get_n_flavor(); ++flavor){
        if(get_flavor_order(flavor) == 0) continue;
        if(c_times[flavor][0].branch != IMAG || a_times[flavor][0].branch != IMAG) return true;
    }
    return false;
}

void Configuration::maximum_real_time(operator_type &op, int &flavor, int &pos) const {
    auto compare_real_time = [](const contour_time& ct_1, const contour_time &ct_2)
    { return (ct_1.branch != IMAG ? ct_1.t : 0) < (ct_2.branch != IMAG ? ct_2.t : 0); };
    contour_time ct_max = {IMAG, 0};
    for(int a = 0; a < get_n_flavor(); ++a){
        auto iter_c = std::max_element(c_times[a].begin(), c_times[a].end(), compare_real_time);
        auto iter_a = std::max_element(a_times[a].begin(), a_times[a].end(), compare_real_time);;
        if(compare_real_time(*iter_a, *iter_c) && compare_real_time(ct_max, *iter_c)){
            ct_max = *iter_c;
            op = creation;
            flavor = a;
            pos = std::distance(c_times[a].begin(), iter_c);
        }
        else if (compare_real_time(*iter_c, *iter_a) && compare_real_time(ct_max, *iter_a)){
            ct_max = *iter_a;
            op = annihilation;
            flavor = a;
            pos = std::distance(a_times[a].begin(), iter_a);
        }
    }
}

std::ostream &operator<<(std::ostream &os, const Configuration &configuration) {
    using std::endl;
    os << "*** Configuration ***" << endl;
    os << "k: " << configuration.k << endl;
    for(int f = 0; f < configuration.get_n_flavor(); ++f){
        os << " flavor no. " << f << ": " << configuration.get_flavor_order(f) << endl;
    }
    os << "abs_weight: " << configuration.abs_weight << endl;
    os << "external_phase: " << std::arg(configuration.external_sign) / PI << " pi" << endl;
    os << "internal_phase: " << std::arg(configuration.internal_sign) / PI << " pi" << endl;
    os << "Initial n: " << endl;
    for(int a = 0; a < configuration.get_n_flavor(); ++a){
        os << " flavor no. " << a << ": " << configuration.get_initial_n(a);
        os << endl;
    }
    if(configuration.k > 0){
        os << "Annihilation times: " << endl;
        int n_flavor = configuration.get_n_flavor();
        for(int a = 0; a < n_flavor; ++a){
            os << " flavor no. " << a << ": ";
            for(const auto &t : configuration.a_times[a]){
                os << t.s(configuration.t_max) << " ";
            }
            os << endl;
        }
        os << "Creation times: " << endl;
        for(int a = 0; a < n_flavor; ++a){
            os << " flavor no. " << a << ": ";
            for(const auto &t : configuration.c_times[a]){
                os << t.s(configuration.t_max) << " ";
            }
            os << endl;
        }
        os << "Segments: " << endl;
        for(int a = 0; a < n_flavor; ++a){
            auto segments = configuration.get_segments(a);
            os << " flavor no. " << a << ": ";
            for(const auto &seg : segments){
                os << "(" << seg.start_time().s(configuration.t_max) << ", "
                        << seg.end_time().s(configuration.t_max) << ", n = " << seg.n() << ") ";
            }
            os << endl;

        }
        os << "(t_max: " << configuration.t_max << ", beta: " << configuration.beta << ")" << endl;
    }
    return os;
}



