#pragma once

#include "./enums/population_type.h"
#include "./structs/member.h"

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void initialize_population(PopulationType population_type, std::uint32_t size);

private:
    Model& model_ref;
    std::vector<Member> population;

    void random_initialization();
    void neighbour_initialization();
};
