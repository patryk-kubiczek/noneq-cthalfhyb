#pragma once

#include "../Moves.h"

template <unsigned N>
class InsertSegment : public Moves::MoveBase {
public:
    using MoveBase::MoveBase;
    std::string name() const override { return (N == 0 ? "Segment insertion " : "Antisegment insertion"); }
    bool execute(int flavor) override;
private:
    contour_time t_start, t_end;
    double t_end_freedom_length;
    int c_index, a_index;
    bool is_not_wound;
    cx_long_double trace_ratio, det_ratio;

    bool randomly_choose_times(int flavor){
        // Select randomly a new start time of a segment(antisegment)
        t_start = rng()->random_contour_time({PLUS, 0}, {IMAG, beta()}, t_max(), beta());
        DEBUG( std::cout << (N == 0 ? "Proposed t_start of a segment: " : "Proposed t_start of an antisegment: ")
                         << t_start.s(t_max()) << std::endl; )

        int k_flavor = conf()->get_flavor_order(flavor);
        if(k_flavor == 0){
            t_end = rng()->random_contour_time({PLUS, 0}, {IMAG, beta()}, t_max(), beta());
            t_end_freedom_length =  beta() + 2 * t_max();
            return true;
        }
        else{
            auto segments = conf()->get_segments(flavor);
            for(const auto &segment : segments){
                if(segment.contains(t_start)){
                    DEBUG( std::cout << (segment.n() == 0 ? "t_start lies in segment: (" : "t_start lies in antisegment: (")
                                     << segment.start_time().s(t_max()) << ", " << segment.end_time().s(t_max())
                                     << ")" << std::endl; )
                    // Reject if start time of a new segment(antisegment) in an already existing segment(antisegment)
                    if(segment.n() == N) return false;
                    else {
                        const contour_time &t_end_max = segment.end_time();
                        t_end = rng()->random_contour_time(t_start, t_end_max, t_max(), beta());
                        t_end_freedom_length  = contour_time::distance(t_end_max, t_start, t_max(), beta());
                        DEBUG( std::cout << "Proposed pair: (t_start, t_end) = (" << t_start.s(t_max()) << ", " << t_end.s(t_max()) << ")" << std::endl;
                               std::cout << "Freedom interval length: " << t_end_freedom_length << std::endl; )
                        return true;
                    }
                }
            }
        }
    }

    bool try_move(int flavor){
        int old_k_flavor = conf()->get_flavor_order(flavor) - 1;
        long double acceptance = (beta() + 2 * t_max()) * t_end_freedom_length / (old_k_flavor + 1)
                                 * std::abs(trace_ratio * det_ratio);
        acceptance = std::min(acceptance, 1.0l);
        return rng()->random_number(0, 1.0) <= acceptance;
    }

    void update_configuration(int flavor, const contour_time &t_c, const contour_time &t_a){
        conf()->increment_k(flavor);
        c_index = get_insertion_position(conf()->get_c_times(flavor), t_c);
        insert_to_vector(conf()->get_c_times(flavor), c_index, t_c);
        a_index = get_insertion_position(conf()->get_a_times(flavor), t_a);
        insert_to_vector(conf()->get_a_times(flavor), a_index, t_a);
    }

    void cancel_update_configuration(int flavor){
        erase_from_vector(conf()->get_a_times(flavor), a_index);
        erase_from_vector(conf()->get_c_times(flavor), c_index);
        conf()->decrement_k(flavor);
    }

#ifdef CTHYB_QMC
    void attempt_insertion(int flavor, const contour_time &t_c, const contour_time &t_a) {
        local_weight()->attempt_add_operator(creation, flavor, t_c);
        local_weight()->attempt_add_operator(annihilation, flavor, t_a);
    }

    void confirm_insertion(int flavor) {
        local_weight()->confirm_add_operator(creation, flavor, c_index);
        local_weight()->confirm_add_operator(annihilation, flavor, a_index);
    }

    void reset_local_weight() {
        local_weight()->cancel_operations();
    }
#endif

#ifdef CTHALFHYB_QMC
    void attempt_insertion(int flavor, const contour_time &t_c, const contour_time &t_a) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_not_wound) lw->attempt_add_pair(N, t_end, t_start);
        else lw->attempt_add_wound_pair(N, t_end, t_start);
    }

    void confirm_insertion(int flavor) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_not_wound) lw->confirm_add_pair(N, c_index, a_index);
        else lw->confirm_add_wound_pair(N);
    }

    void reset_local_weight() {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        lw->cancel_operations();
    }
#endif
};


template<unsigned int N>
bool InsertSegment<N>::execute(int flavor) {
    // Reject the move if the whole contour is one segment (antisegment)
    if(conf()->get_flavor_order(flavor) == 0 && conf()->get_initial_n(flavor) == N) return false;

    // Choose random t_start and t_end from (t_start, t_end_max)
    if(!randomly_choose_times(flavor)) return false;

    // Calculate the ratio of the new and old bath determinant
    const auto& t_c = (N == 0 ? t_end : t_start);
    const auto& t_a = (N == 0 ? t_start : t_end);
    bath_weight()->attempt_insertion(flavor, t_c, t_a);

    det_ratio = bath_weight()->get_det_ratio();

    // Update Configuration
    is_not_wound = t_end > t_start;
    update_configuration(flavor, t_c, t_a);

    // Calculate new_local trace
    bool is_local_weight = local_weight()->type == NUMERIC_LW || local_weight()->type == ANALYTIC_LW;
    if(is_local_weight) attempt_insertion(flavor, t_c, t_a);
    local_weight()->calculate_new_trace();
    trace_ratio = local_weight()->get_trace_ratio();

    // Accept or reject the move
    if(!try_move(flavor)) {
        cancel_update_configuration(flavor);
        if(is_local_weight) reset_local_weight();
        return false;
    }

    // Update Local Weight
    if(is_local_weight) confirm_insertion(flavor);
    local_weight()->finalize_update();

    // Update Bath Weight
    bath_weight()->confirm_insertion(conf()->get_absolute_index(creation, flavor, c_index),
                                     conf()->get_absolute_index(annihilation, flavor, a_index));

#ifdef CTHYB_QMC
    // Dirty fix for instability arising for large beta (2/2)
    if(beta() > 1)
        bath_weight()->generate_from_conf();
#endif

    // Update weight
    conf()->calculate_internal_sign();
    conf()->set_weight(local_weight()->get_w_loc() * cx_long_double{bath_weight()->get_w_bath()});
    return true;
}











