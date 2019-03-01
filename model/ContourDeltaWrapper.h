#pragma once


#ifdef CTHYB_QMC
#include "../model/cthyb/Model.h"

using contour_delta_t = ContourDelta;

#endif

#ifdef CTHALFHYB_QMC
#include "../model/cthalfhyb/Model.h"

using contour_delta_t = ContourDeltaSingle;

#endif

enum bath_weight_type { NUMERIC_BW, ANALYTIC_BW };

class ContourDeltaWrapper {
public:
    static ContourDeltaWrapper create(const Model& model, bath_weight_type type){
#ifdef CTHYB_QMC
        if(type == NUMERIC_BW)
            return {model.delta};
        else // if(type == ANALYTIC_BW)
            return {model.delta_up_analytic, model.delta_down_analytic};
#endif
#ifdef CTHALFHYB_QMC
        if(type == NUMERIC_BW)
            return {model.delta_down};
        else  // if(type == ANALYTIC_BW)
            return {model.delta_down_analytic, model.delta_down_analytic};
#endif
    }
    cx_double operator()(int c_flavor, int a_flavor, const contour_time& ct_1, const contour_time& ct_2) const {
        if(!analytic)
            return (*numeric_delta)(c_flavor, a_flavor, ct_1, ct_2);
        else{
            if(c_flavor == a_flavor){
                if(c_flavor == 0) return (*analytic_delta_up)(ct_1, ct_2);
                if(c_flavor == 1) return (*analytic_delta_down)(ct_1, ct_2);
            }
            else return 0.;
        }
    }

private:
    ContourDeltaWrapper(const contour_delta_t &delta)
            : numeric_delta(&delta), analytic_delta_down(nullptr), analytic_delta_up(nullptr), analytic(false) {}
    ContourDeltaWrapper(const AnalyticDelta &delta_up, const AnalyticDelta &delta_down)
            : numeric_delta(nullptr), analytic_delta_up(&delta_up), analytic_delta_down(&delta_down), analytic(true) {}
    const contour_delta_t *numeric_delta;
    const AnalyticDelta *analytic_delta_down;
    const AnalyticDelta *analytic_delta_up;
    bool analytic;
};
