#pragma once

#include "./enums/population_type.h"
#include "./structs/member.h"

#include <tuple>
#include <vector>

using OffspringsPair = std::pair<std::vector<Node>, std::vector<Node>>;

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void initialize_population(PopulationType population_type, std::size_t size);
    void print_population();

    std::vector<Member> tournament_selection(std::size_t subgroup_size);
    OffspringsPair ordered_crossover(const Member& parent1, const Member& parent2);
    OffspringsPair ordered_crossover();

// private:
    Model& model_ref;
    std::vector<Member> population;

    void random_initialization();
    void neighbour_initialization();
    void mixed_initialization();
};
