#pragma once

#ifdef CTHALFHYB_QMC
#include "cthalfhyb/LocalWeight.h"
#include "../discrete_bath/DiscreteBathWeight.h"

#ifndef LOCAL_WEIGHT_EVALUATION
#define LOCAL_WEIGHT_EVALUATION LocalWeightBase
#endif

using local_weight_t = LOCAL_WEIGHT_EVALUATION;

#endif

#ifdef CTHYB_QMC
#include "cthyb/LocalWeight.h"

using local_weight_t = LocalWeight;

#endif

template <typename LW = LocalWeightBase>
class LocalWeightWrapper {
public:
    LocalWeightWrapper(const Model &m, const Configuration &conf, local_weight_type type)
            : local_weight_ptr(create_ptr(m, conf, type)) {}
    const LW& get() const { return *local_weight_ptr; }
    LW& get() { return *local_weight_ptr; }
private:
    std::unique_ptr<LW> local_weight_ptr;

    template<bool cond>
    using unique_ptr_LW_enable_if = typename std::enable_if_t<cond, std::unique_ptr<LW>>;

    template <typename LWeight = LW>
    static unique_ptr_LW_enable_if<std::is_same<LWeight, LocalWeightBase>::value>
            create_ptr(const Model &m, const Configuration &conf, local_weight_type type) {
        switch(type) {
            case NUMERIC_LW:
                return std::make_unique<LocalWeight>(m, conf, type);
            case ANALYTIC_LW:
                return std::make_unique<LocalWeight>(m, conf, type);
#ifdef CTHALFHYB_QMC
            case DISCRETE_BATH:
            return std::make_unique<DiscreteBathWeight>(m, conf, type);
#endif
        }
    }

    template <typename LWeight = LW>
    static unique_ptr_LW_enable_if<!std::is_same<LWeight, LocalWeightBase>::value>
            create_ptr(const Model &m, const Configuration &conf, local_weight_type type) {
        return std::make_unique<LW>(m, conf, type);
    }
};

