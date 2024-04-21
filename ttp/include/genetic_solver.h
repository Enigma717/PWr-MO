#pragma once

#include "./enums/population_type.h"
#include "./structs/solution.h"

#include <tuple>
#include <vector>
#include <fstream>

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void print_population();
    void evaluate_population();
    const Solution& solve();

private:
    Model& model_ref;
    std::vector<Solution> population;
    Solution best_solution;
    Solution worst_solution;
    double avg_fitness;

    std::size_t fitness_evaluations {0};
    std::size_t generation_number {0};
    std::size_t population_size {200};
    std::size_t subgroup_size {7};
    double crossing_probability {0.8};
    double mutation_probability {0.1};

    void initialize_population(PopulationType population_type, std::size_t size);
    void random_initialization();
    void neighbour_initialization();
    void mixed_initialization();

    std::vector<Solution> tournament_selection(std::size_t subgroup_size);
    std::vector<Solution> process_crossover(const std::vector<Solution>& tournament_winners);
    void evolve_population(
        const std::vector<Solution>& parents, const std::vector<Solution>& offsprings);
    void process_mutation();

    Solution order_crossover(
        const Solution& first_parent,
        const Solution& second_parent,
        const std::size_t first_crossing_point,
        const std::size_t second_crossing_point);
    Solution create_new_solution(const std::vector<Node>& route);
    double fitness_evaluation(Solution& solution);
};
