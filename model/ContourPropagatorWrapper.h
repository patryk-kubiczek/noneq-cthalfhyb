#pragma once

#include "../local_weight/LocalWeightBase.h" // to get definition of local_weight_type

#ifdef CTHYB_QMC
#include "../model/cthyb/Model.h"

class ContourPropagatorWrapper{
public:
    static ContourPropagatorWrapper create(const Model& m, local_weight_type type){
        if(type == NUMERIC_LW)
            return {m.bare_propagator};
        else // if(type == ANALYTIC_LW)
            return {m.bare_propagator_analytic};
    }
    arma::cx_mat operator()(const contour_time& ct_1, const contour_time& ct_2) const {
        return (analytic ? (*analytic_p)(ct_1, ct_2) : (*numeric_p)(ct_1, ct_2));
    }
private:
    ContourPropagatorWrapper(const ContourPropagator &p)
            : numeric_p(&p), analytic_p(nullptr), analytic(false) {}
    ContourPropagatorWrapper(const AnalyticP &p)
            : numeric_p(nullptr), analytic_p(&p), analytic(true) {}
    const ContourPropagator *numeric_p;
    const AnalyticP *analytic_p;
    bool analytic;
};

#endif

#ifdef CTHALFHYB_QMC
#include "../model/cthalfhyb/Model.h"

template<unsigned n>
class ContourUWrapper{
    template<bool cond>
    using ContourUWrapper_enable_if = typename std::enable_if_t<cond, ContourUWrapper>;
public:
    template<unsigned n_down = n>
    static ContourUWrapper_enable_if<n_down == 0> create(const Model& m, local_weight_type type){
        if(type == NUMERIC_LW)
            return {m.u_0};
        else // if(type == ANALYTIC_LW)
            return {m.u_0_analytic};
    }
    template<unsigned n_down = n>
    static ContourUWrapper_enable_if<n_down == 1> create(const Model& m, local_weight_type type){
        if(type == NUMERIC_LW)
            return {m.u_1};
        else // if(type == ANALYTIC_LW)
            return {m.u_1_analytic};
    }
    arma::cx_mat operator()(const contour_time& ct_1, const contour_time& ct_2) const {
        return (analytic ? (*analytic_u)(ct_1, ct_2) : (*numeric_u)(ct_1, ct_2));
    }

private:
    ContourUWrapper(const ContourPropagator &u)
            : numeric_u(&u), analytic_u(nullptr), analytic(false) {}
    ContourUWrapper(const AnalyticU<n> &u)
            : numeric_u(nullptr), analytic_u(&u), analytic(true) {}
    const ContourPropagator *numeric_u;
    const AnalyticU<n> *analytic_u;
    bool analytic;
};

class ContourPhi1Wrapper{
public:
    static ContourPhi1Wrapper create(const Model& m, local_weight_type type){
        if(type == NUMERIC_LW)
            return {m.phi_1};
        else // if(type == ANALYTIC_LW)
            return {m.phi_1_analytic};
    }
    cx_double operator()(const contour_time& ct_1, const contour_time& ct_2) const {
        return (analytic ? (*analytic_phi)(ct_1, ct_2) : (*numeric_phi)(ct_1, ct_2));
    }
private:
    ContourPhi1Wrapper(const ContourPhi &phi)
            : numeric_phi(&phi), analytic_phi(nullptr), analytic(false) {}
    ContourPhi1Wrapper(const AnalyticPhi1 &phi)
            : numeric_phi(nullptr), analytic_phi(&phi), analytic(true) {}
    const ContourPhi *numeric_phi;
    const AnalyticPhi1 *analytic_phi;
    bool analytic;
};

#endif




