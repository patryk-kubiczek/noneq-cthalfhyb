#pragma once

#include "../Moves.h"

template <unsigned N>
class RemoveSegment : public Moves::MoveBase {
public:
    using MoveBase::MoveBase;
    std::string name() const override { return (N == 0 ? "Segment removal" : "Antisegment removal"); }
    bool execute(int flavor) override;
private:
    Configuration::Segment segment;
    double t_end_freedom_length;
    bool is_not_wound;
    cx_long_double trace_ratio, det_ratio;

    void randomly_choose_segment(int flavor){
        auto segments = conf()->get_segments(flavor);
        int i = rng()->random_int(0, segments.size() / 2 - 1);
        if(segments[0].n() == N) i = 2 * i;
        else i = 2 * i + 1;
        segment = segments[i];
        assert(segment.n() == N);
        t_end_freedom_length = contour_time::distance(segments[(i + 1) % segments.size()].end_time(),
                                       segment.start_time(), t_max(), beta());
        DEBUG( std::cout << "Trying to remove " << (N == 0 ? "segment (" : "antisegment (")
                << segment.start_time().s(t_max()) << ", " << segment.end_time().s(t_max()) << ")" << std::endl; )
        
    }

    bool try_move(int flavor){
        int old_k_flavor = conf()->get_flavor_order(flavor) + 1;
        long double acceptance = old_k_flavor / ((beta() + 2 * t_max()) * t_end_freedom_length)
                                 * std::abs(trace_ratio * det_ratio);
        acceptance = std::min(acceptance, 1.0l);
        return rng()->random_number(0, 1.0) <= acceptance;
    }

    void update_configuration(int flavor){
        int initial_n = conf()->get_initial_n(flavor);
        erase_from_vector(conf()->get_c_times(flavor), segment.c_position());
        erase_from_vector(conf()->get_a_times(flavor), segment.a_position());
        conf()->decrement_k(flavor);
        if(is_not_wound) conf()->set_initial_n(flavor, initial_n);
        else conf()->set_initial_n(flavor, (initial_n + 1) % 2);
    }

    void cancel_update_configuration(int flavor){
        conf()->increment_k(flavor);
        insert_to_vector(conf()->get_c_times(flavor), segment.c_position(), segment.t_c());
        insert_to_vector(conf()->get_a_times(flavor), segment.a_position(), segment.t_a());
    }

#ifdef CTHYB_QMC
    void attempt_removal(int flavor) {
        local_weight()->attempt_remove_operator(creation, flavor, segment.c_position());
        local_weight()->attempt_remove_operator(annihilation, flavor, segment.a_position());
    }

    void confirm_removal(int flavor) {
        local_weight()->confirm_remove_operator(creation, flavor, segment.c_position());
        local_weight()->confirm_remove_operator(annihilation, flavor, segment.a_position());
    }

    void reset_local_weight() {
        local_weight()->cancel_operations();
    }
#endif

#ifdef CTHALFHYB_QMC
    void attempt_removal(int flavor) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_not_wound) lw->attempt_remove_pair(N, (N == 0 ? segment.c_position() : segment.a_position()));
        else lw->attempt_remove_wound_pair(N);
    }

    void confirm_removal(int flavor) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_not_wound)lw->confirm_remove_pair(N, segment.c_position(), segment.a_position());
        else lw->confirm_remove_wound_pair(N);
    }

    void reset_local_weight() {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        lw->cancel_operations();
    }
#endif
};

template<unsigned int N>
bool RemoveSegment<N>::execute(int flavor) {
    // Reject the move if the flavor order = 0
    int k_flavor = conf()->get_flavor_order(flavor);
    if(k_flavor == 0) return false;

    // Choose random c_index and a corresponding a_index
    randomly_choose_segment(flavor);

    // Calculate the ratio of the new and old bath determinant
    int i = conf()->get_absolute_index(creation, flavor, segment.c_position());
    int j = conf()->get_absolute_index(annihilation, flavor, segment.a_position());
    bath_weight()->attempt_removal(i, j);
    det_ratio = bath_weight()->get_det_ratio();

    // Update Configuration
    is_not_wound = !segment.is_wound();
    update_configuration(flavor);

    // Calculate the new local trace
    bool is_local_weight = local_weight()->type == NUMERIC_LW || local_weight()->type == ANALYTIC_LW;
    if(is_local_weight) attempt_removal(flavor);
    local_weight()->calculate_new_trace();
    trace_ratio = local_weight()->get_trace_ratio();

    // Accept or reject the move
    if(!try_move(flavor)) {
        cancel_update_configuration(flavor);
        if(is_local_weight) reset_local_weight();
        return false;
    }

    // Update Local Weight
    if(is_local_weight) confirm_removal(flavor);
    local_weight()->finalize_update();

    // Update Bath Weight
    bath_weight()->confirm_removal(i, j);

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
