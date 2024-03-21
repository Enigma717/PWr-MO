#pragma once

#include "./enums/population_type.h"
#include "./structs/human.h"

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void initialize_population(PopulationType population_type, std::uint32_t size);

private:
    Model& model_ref;
    std::vector<Human> population;

    void random_initialization();
    void neighbour_initialization();
};
