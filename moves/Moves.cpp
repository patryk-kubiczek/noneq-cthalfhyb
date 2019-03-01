#include "Moves.h"

#include "moves/ToggleFullEmptyLine.h"
#include "moves/InsertSegment.h"
#include "moves/RemoveSegment.h"
#include "moves/ShiftSegmentEnd.h"

//#include "old_moves/InsertSegment.h"
//#include "old_moves/InsertAntisegment.h"
//#include "old_moves/RemoveSegment.h"
//#include "old_moves/RemoveAntisegment.h"


Moves::Moves(const Model& model, Configuration &conf, local_weight_t &local_weight,
             BathWeight &bath_weight, RandomNumberGenerator &rng)
        : t_max(model.t_max), beta(model.beta), conf(&conf), local_weight(&local_weight),
          bath_weight(&bath_weight), rng(&rng), model(&model) {

    initialize_moves();
    sweep_size = 0;
    for(const auto& m : moves)
        sweep_size += m->multiplicity;
}

void Moves::initialize_moves() {
    using std::make_unique;
    int index = 0;

    moves.emplace_back(make_unique<ToggleFullEmptyLine>(*this, 1, index++));
    moves.emplace_back(make_unique<InsertSegment<0>>(*this, 1, index++));
    moves.emplace_back(make_unique<InsertSegment<1>>(*this, 1, index++));
    moves.emplace_back(make_unique<RemoveSegment<0>>(*this, 1, index++));
    moves.emplace_back(make_unique<RemoveSegment<1>>(*this, 1, index++));
    moves.emplace_back(make_unique<ShiftSegmentEnd>(*this, 4, index));
}

Moves::MoveBase& Moves::random_move() const {
    int random_i = rng->random_int(0, sweep_size - 1);
    int i = 0;
    for(const auto& m : moves){
        for(int j = 0; j < m->multiplicity; ++j){
            if(i == random_i) return (*m);
            else ++i;
        }
    }
}

std::ostream &operator<<(std::ostream &os, const Moves &moves) {
    using std::endl;
    os << "*** Moves ***" << endl;
    os << "Sweep:" << endl;
    for(const auto &m : moves.moves){
        os << " " << m->multiplicity << " x [" << m->name() << "]" << endl;
    }
    os << "Sweep size: " << moves.sweep_size << endl;
    return os;
}




