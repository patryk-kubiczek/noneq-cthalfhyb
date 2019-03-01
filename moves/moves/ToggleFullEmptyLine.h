#pragma once

#include "../Moves.h"

class ToggleFullEmptyLine : public Moves::MoveBase {
public:
    using MoveBase::MoveBase;
    std::string name() const override { return "Toggle full/empty line"; }
    bool execute(int flavor) override;

private:
#ifdef CTHALFHYB_QMC
    void attempt_switch_n(int new_n_down) {
        auto* lw = static_cast<LocalWeight*>(local_weight());
        lw->attempt_switch_n(conf()->get_initial_n());
    }

    void reset_local_weight(){
        auto* lw = static_cast<LocalWeight*>(local_weight());
        lw->cancel_operations();
    }
#endif
#ifdef CTHYB_QMC
    void attempt_switch_n(int new_n_down) {}
    void reset_local_weight() {}
#endif
};

bool ToggleFullEmptyLine::execute(int flavor) {
    if(conf()->get_flavor_order(flavor) != 0)
        return false;

    int old_initial_n = conf()->get_initial_n(flavor);
    conf()->set_initial_n(flavor, (old_initial_n + 1) % 2);

    bool is_local_weight = local_weight()->type == NUMERIC_LW || local_weight()->type == ANALYTIC_LW;
    if(is_local_weight) attempt_switch_n(conf()->get_initial_n());
    local_weight()->calculate_new_trace();

    auto acceptance = std::min(std::abs(local_weight()->get_trace_ratio()), 1.0l);
    if(rng()->random_number(0, 1.0) > acceptance) {
        conf()->set_initial_n(flavor, old_initial_n);
        if(is_local_weight) reset_local_weight();
        return false;
    }

    local_weight()->finalize_update();
    conf()->set_weight(local_weight()->get_w_loc() * cx_long_double{bath_weight()->get_w_bath()});

    return true;
}

