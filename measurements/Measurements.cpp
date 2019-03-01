#include "Measurements.h"

#include "measurements/MeanComplexSign.h"
#include "measurements/MeanOrder.h"
#include "measurements/OrderHistogram.h"
#include "measurements/InitialOccupancies.h"
#include "measurements/Occupancies.h"
#include "measurements/GreenFunctions.h"
#include "measurements/EquilibriumGreenFunctions.h"
#include "measurements/MatsubaraGreenFunctions.h"
#include "measurements/TimesStatistics.h"
#include "measurements/DoubleOccupancy.h"

#ifdef CTHALFHYB_QMC
#include "measurements/SpinUp.h"
#endif


void Measurements::initialize_measurements() {
    acceptance_rates.resize(moves->number_of_moves(), 0);

    using std::make_unique;

    // ********* First measurement has to be always MeanComplexSign! **********
    measurements.emplace_back(make_unique<MeanComplexSign>(*this));
    // ************************************************************************
    measurements.emplace_back(make_unique<MeanOrder>(*this));
    measurements.emplace_back(make_unique<OrderHistogram>(*this));
    measurements.emplace_back(make_unique<InitialOccupancies>(*this));
    if(t_max > 0)
        measurements.emplace_back(make_unique<Occupancies>(*this));
#ifdef CTHALFHYB_QMC
    measurements.emplace_back(make_unique<SpinUp>(*this));
#endif
#ifdef CTHYB_QMC
    measurements.emplace_back(make_unique<DoubleOccupancy>(*this));
#endif
    measurements.emplace_back(make_unique<GreenFunctions>(*this));
    //measurements.emplace_back(make_unique<EquilibriumGreenFunctions>(*this));
    if(n_imag_times > 1)
        measurements.emplace_back(make_unique<MatsubaraGreenFunctions>(*this));
    measurements.emplace_back(make_unique<TimeStatistics>(*this));
    for(auto &m : measurements){
        m->reset();
    }
}


void Measurements::increment_exception_counter(const std::string &exception_name) {
    if(exception_counter.count(exception_name) == 0)
        exception_counter[exception_name] = 1;
    else
        ++exception_counter[exception_name];
}


void Measurements::do_measurements() {
    if(step_counter % block_size != 0) return;
    ++measurement_counter;
    for(auto &m : measurements){
        m->measure();
    }
}


void Measurements::reset_measurements() {
    measurement_counter = 0;
    step_counter = 0;
    exception_counter.clear();
    for(auto &acc : acceptance_rates){
        acc = 0;
    }
    for(auto &m : measurements){
        m->reset();
    }
}


cx_double Measurements::mean_sign() const{
    auto &m = static_cast<MeanComplexSign&>(*(measurements.front()));
    return m.result();
}


double Measurements::mean_acceptance() const{
    long long result = 0;
    for(int i = 0; i < moves->number_of_moves(); ++i) {
        result += acceptance_rates[i];
    }
    return double(result) / step_counter;
}


void Measurements::print_acceptance_rates(std::ostream &os) const{
    using std::endl;
    os << "Acceptancies:" << endl;
    for(int i = 0; i < moves->number_of_moves(); ++i){
        os << " [" << moves->move(i).name() << "]: "
           << double(acceptance_rates[i]) / (double(step_counter) * moves->move(i).multiplicity / moves->get_sweep_size()) << endl;
    }
}

std::ostream &operator<<(std::ostream &os, const Measurements &meas) {
    using std::endl;
    os << "*** Measurements ***" << endl;
    os << "MC steps: " << meas.step_counter << endl;
    os << "MC measurements: " << meas.measurement_counter << endl;
    meas.print_acceptance_rates(os);
    os << "Mean acceptance: " << meas.mean_acceptance() << endl;
    if(!meas.exception_counter.empty()){
        os << "Exceptions:" << endl;
        for(const auto& el : meas.exception_counter){
            os << " [" << el.first << "]: " << el.second << endl;
        }
    }
    for(const auto& m : meas.measurements) {
        m->print_normalized(os);
    }
    return os;
}

