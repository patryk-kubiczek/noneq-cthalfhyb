#pragma once

#include <ostream>
#include <type_traits>

#include "../configuration/Configuration.h"
#include "../bath_weight/BathWeight.h"
#include "../random_number_generator/RandomNumberGenerator.h"
#include "../local_weight/LocalWeightWrapper.h"

class Moves {
public:
    Moves(const Model &model, Configuration &conf, local_weight_t &local_weight, BathWeight &bath_weight, 
          RandomNumberGenerator &rng);

    class MoveBase{
    public:
        MoveBase(Moves &moves, int multiplicity, int index) : moves(&moves), multiplicity(multiplicity), index(index) {}
        virtual bool execute(int flavor) = 0;
        virtual std::string name() const = 0;
        const int multiplicity;
        const int index;
        virtual ~MoveBase() = default;
    protected:
        Configuration* conf() { return moves->conf; }
        local_weight_t* local_weight() { return moves->local_weight; }
        BathWeight* bath_weight() { return moves->bath_weight; }
        RandomNumberGenerator* rng() { return moves->rng; }
        // ***
        const Model* model() { return moves->model; }
        // ***
        double t_max() { return moves->t_max; }
        double beta() { return moves->beta; }
    private:
        Moves *moves;
    };

    MoveBase& move(int i) const { return *(moves[i]); }
    MoveBase& random_move() const;
    int number_of_moves() const { return moves.size(); }
    int get_sweep_size() const { return sweep_size; }

    friend std::ostream &operator<<(std::ostream &os, const Moves &moves);

private:
    double t_max;
    double beta;

    // ***
    const Model *model;
    // ***
    Configuration *conf;
    local_weight_t *local_weight;
    BathWeight *bath_weight;
    RandomNumberGenerator *rng;

    std::vector<std::unique_ptr<MoveBase>> moves;
    int sweep_size;
    
    void initialize_moves();
};