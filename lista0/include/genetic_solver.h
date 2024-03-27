#pragma once

#include "./enums/population_type.h"
#include "./structs/member.h"

#include <tuple>
#include <vector>

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void print_population();
    void evaluate_population();
    void solve();

    void initialize_population(PopulationType population_type, std::size_t size);
    std::vector<Member> tournament_selection(std::size_t subgroup_size);
    std::vector<Member> process_crossover(const std::vector<Member>& tournament_winners);
    void evolve_population(
        const std::vector<Member>& parents, const std::vector<Member>& offsprings);
    void process_mutation();

private:
    Model& model_ref;
    std::vector<Member> population;
    Member best_member;

    std::size_t generation_number {0};
    std::size_t population_size {200};
    std::size_t subgroup_size {7};
    double crossing_probability {0.8};
    double mutation_probability {0.1};

    void random_initialization();
    void neighbour_initialization();
    void mixed_initialization();

    Member order_crossover(
        const Member& first_parent,
        const Member& second_parent,
        const std::size_t first_crossing_point,
        const std::size_t second_crossing_point);
};
