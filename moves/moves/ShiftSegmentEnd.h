#pragma once

#include "../Moves.h"

class ShiftSegmentEnd : public Moves::MoveBase {
public:
    using MoveBase::MoveBase;
    std::string name() const override { return "Shift segment/antisegment end"; }
    bool execute(int flavor) override;
private:
    int segment_position;
    Configuration::Segment segment;
    contour_time new_ct, old_ct;
    bool is_old_not_wound, is_new_not_wound;
    operator_type op_type;
    int op_index, new_op_index;
    int winding_direction;
    cx_long_double trace_ratio, det_ratio;

    void randomly_choose_segment_and_new_time(int flavor) {
        auto segments = conf()->get_segments(flavor);
        segment_position = rng()->random_int(0, segments.size() - 1);
        segment = segments[segment_position];
        old_ct = segment.end_time();
        DEBUG( std::cout << "Limits for new end time: [" << segment.start_time().s(t_max()) << ", "
                         << segments[(segment_position + 1) % segments.size()].end_time().s(t_max()) <<  "]" << std::endl; )

        if(conf()-> get_flavor_order(flavor) > 1)
            new_ct = rng()->random_contour_time(segment.start_time(), segments[(segment_position + 1) % segments.size()].end_time(),
                                           t_max(), beta());
        else
            new_ct = rng()->random_contour_time({PLUS, 0}, {IMAG, beta()}, t_max(), beta());
        DEBUG( std::cout << "Old end time: " << old_ct.s(t_max()) << ", new end time: " << new_ct.s(t_max()) << std::endl; )
        is_old_not_wound = !segment.is_wound();
        is_new_not_wound = segment.start_time() < new_ct;
        op_index = (segment.n() == 0 ? segment.c_position() : segment.a_position());
        op_type = (segment.n() == 0 ? creation : annihilation);
    }

    bool try_move(int flavor) {
        long double acceptance = std::abs(trace_ratio * det_ratio);
        acceptance = std::min(acceptance, 1.0l);
        return rng()->random_number(0, 1.0) <= acceptance;
    }

    void find_new_op_index(int flavor) {
        int last_position = 2 * conf()->get_flavor_order(flavor) - 1;
        if(is_new_not_wound == is_old_not_wound){
            new_op_index = op_index;
            winding_direction = 0;
        }
        else if(segment_position == last_position - 1){
            new_op_index = 0;
            winding_direction = 1;
        }
        else if(segment_position == last_position){
            new_op_index = conf()->get_flavor_order(flavor) - 1;
            winding_direction = -1;
        }
        DEBUG( else{ std::cout << "I shouldn't be here. " << std::endl;
                     std::cout << "Last position: " << last_position << ", seg. pos: " << segment_position << std::endl; } )
        DEBUG( std::cout << "Winding direction: " << winding_direction << std::endl; )
    }

    void update_configuration(int flavor){
        if(new_op_index == op_index){
            if(segment.n() == 0){
                conf()->get_contour_time(creation, flavor, op_index) = new_ct;
            }
            else{
                conf()->get_contour_time(annihilation, flavor, op_index) = new_ct;
            }
        }
        else{
            if(segment.n() == 0){
                erase_from_vector(conf()->get_c_times(flavor), op_index);
                insert_to_vector(conf()->get_c_times(flavor), new_op_index, new_ct);
            }
            else{
                erase_from_vector(conf()->get_a_times(flavor), op_index);
                insert_to_vector(conf()->get_a_times(flavor), new_op_index, new_ct);
            }
        }
    }
    void cancel_update_configuration(int flavor){
        if(new_op_index == op_index){
            if(segment.n() == 0){
                conf()->get_contour_time(creation, flavor, op_index) = old_ct;
            }
            else{
                conf()->get_contour_time(annihilation, flavor, op_index) = old_ct;
            }
        }
        else{
            if(segment.n() == 0){
                erase_from_vector(conf()->get_c_times(flavor), new_op_index);
                insert_to_vector(conf()->get_c_times(flavor), op_index, old_ct);
            }
            else{
                erase_from_vector(conf()->get_a_times(flavor), new_op_index);
                insert_to_vector(conf()->get_a_times(flavor), op_index, old_ct);
            }
        }
    }

#ifdef CTHYB_QMC
    void attempt_shift(int flavor) {
        local_weight()->attempt_remove_operator(op_type, flavor, op_index);
        local_weight()->attempt_add_operator(op_type, flavor, new_ct);
    }

    void confirm_shift(int flavor) {
        local_weight()->confirm_remove_operator(op_type, flavor, op_index);
        local_weight()->confirm_add_operator(op_type, flavor, new_op_index);
    }

    void reset_local_weight() {
        local_weight()->cancel_operations();
    }
#endif

#ifdef CTHALFHYB_QMC
    void attempt_shift(int flavor) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_old_not_wound) lw->attempt_remove_pair(segment.n(), op_index);
        else lw->attempt_remove_wound_pair(segment.n());
        if(is_new_not_wound) lw->attempt_add_pair(segment.n(), new_ct, segment.start_time());
        else lw->attempt_add_wound_pair(segment.n(), new_ct, segment.start_time());
    }

    void confirm_shift(int flavor) {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        if(is_old_not_wound) lw->confirm_remove_pair(segment.n(), segment.c_position(), segment.a_position());
        else lw->confirm_remove_wound_pair(segment.n());
        if(is_new_not_wound){
            int new_c_index = (segment.n() == 0 ? new_op_index : segment.c_position());
            int new_a_index = (segment.n() == 0 ? segment.a_position() : new_op_index);
            lw->confirm_add_pair(segment.n(), new_c_index, new_a_index);
        }
        else
            lw->confirm_add_wound_pair(segment.n());
    }

    void reset_local_weight() {
        auto *lw = static_cast<LocalWeight*>(local_weight());
        lw->cancel_operations();
    }
#endif
};

bool ShiftSegmentEnd::execute(int flavor) {
    // Reject the move if the flavor order = 0
    if(conf()->get_flavor_order(flavor) == 0) return false;

    // Choose a random segment
    randomly_choose_segment_and_new_time(flavor);

    // Calculate the ratio of the new and old bath determinant
    if(segment.n() == 0)
        bath_weight()->attempt_shift_t_c(flavor, op_index, new_ct);
    else
        bath_weight()->attempt_shift_t_a(flavor, op_index, new_ct);

    det_ratio = bath_weight()->get_det_ratio();

    // Find new operator index in case of winding
    find_new_op_index(flavor);

    // Update Configuration
    update_configuration(flavor);

    // Calculate the new local trace
    bool is_local_weight = local_weight()->type == NUMERIC_LW || local_weight()->type == ANALYTIC_LW;
    if(is_local_weight) attempt_shift(flavor);
    local_weight()->calculate_new_trace();
    trace_ratio = local_weight()->get_trace_ratio();

    // Accept or reject the move
    if(!try_move(flavor)) {
        cancel_update_configuration(flavor);
        if(is_local_weight) reset_local_weight();
        return false;
    }

    // Update Local Weight
    if(is_local_weight) confirm_shift(flavor);
    local_weight()->finalize_update();

    // Update Bath Weight
    if(segment.n() == 0)
        bath_weight()->confirm_shift_t_c(flavor, op_index, winding_direction);
    else
        bath_weight()->confirm_shift_t_a(flavor, op_index, winding_direction);

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
