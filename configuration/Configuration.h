#pragma once

#include <vector>
#include <armadillo>
#include <ostream>

#include "../random_number_generator/RandomNumberGenerator.h"
#include "../model/contour_time.h"
#include "../auxiliary_functions/constants.h"

enum operator_type{creation, annihilation};

class Configuration {
    using times_container_t = std::vector<contour_time>;
public:
    Configuration(double t_max, double beta, int n_flavor);
    Configuration(const Configuration &conf) = default;
    Configuration& operator=(Configuration &&conf) = default;
    Configuration& operator=(const Configuration &conf) = default;

    // Constant access methods
    int get_k() const { return k; }
    long double get_abs_weight() const { return abs_weight; }
    cx_double get_sign() const { return internal_sign * external_sign; }
    const auto& get_c_times(int flavor = 0) const { return c_times[flavor]; }
    const auto& get_a_times(int flavor = 0) const { return a_times[flavor]; }

    int get_initial_n(int flavor = 0) const;
    void set_initial_n(int flavor, int n) { initial_n[flavor] = n; }

    class Segment {
    public:
        Segment() = default;
        Segment(int c_position, int a_position, const contour_time &t_c, const contour_time &t_a, int n, int flavor)
                : c_position_(c_position), a_position_(a_position), t_c_(t_c), t_a_(t_a), n_(n), flavor(flavor) {}
        const contour_time& t_c() const { return t_c_; }
        const contour_time& t_a() const { return t_a_; }
        const contour_time& start_time() const { return (n_ == 0 ? t_a() : t_c()); }
        const contour_time& end_time() const { return (n_ == 0 ? t_c() : t_a()); }
        const int n() const { return n_; }
        const int c_position() const { return c_position_; }
        const int a_position() const { return a_position_; }
        bool is_wound() const { return start_time() > end_time(); }
        bool contains(const contour_time& ct) const {
            if(!is_wound()) return start_time() < ct && ct < end_time();
            else return start_time() < ct || ct < end_time();
        }
    private:
        int n_;
        int c_position_;
        int a_position_;
        contour_time t_c_;
        contour_time t_a_;
        int flavor;
    };


    std::vector<Segment> get_segments(int flavor = 0) const;
    std::vector<contour_time> get_contour_times(int flavor = 0) const;

    // Absolute index = position among all the operators of a given type
    int get_absolute_index(operator_type op_type, int flavor, int relative_index) const;
    // Global index = position among all the operators
    int get_global_index(operator_type op_type, int flavor, int relative_index) const;
    // Relative index = position among all the operators of given type and flavor
    int get_relative_index(int global_index, operator_type op_type, int &flavor) const;

    const contour_time& get_contour_time(operator_type op_type, int flavor, int relative_index) const {
        return (op_type == annihilation ? a_times : c_times)[flavor][relative_index];
    }
    int get_n_flavor() const{
        assert(c_times.size() == a_times.size());
        return c_times.size();
    }
    int get_flavor_order(int flavor) const{
        return flavor_order[flavor];
    }
    bool has_operators_on_real_branch() const;
    void maximum_real_time(operator_type &op, int &flavor, int &pos) const;

    // Non-const access methods
    auto& get_c_times(int flavor = 0) { return c_times[flavor]; }
    auto& get_a_times(int flavor = 0) { return a_times[flavor]; }
    contour_time& get_contour_time(operator_type op_type, int flavor, int relative_index) {
        return (op_type == annihilation ? a_times : c_times)[flavor][relative_index];
    }
    void increment_k(int flavor) { ++k; ++flavor_order[flavor]; }
    void decrement_k(int flavor) { --k; --flavor_order[flavor]; }
    void set_weight(cx_long_double weight) {
        abs_weight = std::abs(weight);
        external_sign = weight / abs_weight;
    }

    void calculate_internal_sign();
    void randomly_generate(int k, RandomNumberGenerator &rng, bool segment_like = true);

    friend std::ostream &operator<<(std::ostream &os, const Configuration &configuration);

private:
    double t_max;
    double beta;

    int k = 0;
    long double abs_weight = 1;
    cx_double external_sign = 1.;
    cx_double internal_sign = 1.;
    std::vector<times_container_t> c_times;
    std::vector<times_container_t> a_times;
    std::vector<int> initial_n;

    std::vector<int> flavor_order;
    
    void randomly_generate_general(int k_max, RandomNumberGenerator &rng);
    void randomly_generate_segment_like(int k_max, RandomNumberGenerator &rng);
};