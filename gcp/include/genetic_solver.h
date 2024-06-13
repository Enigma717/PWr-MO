#pragma once

#include "./enums/population_type.h"
#include "./structs/solution.h"
#include "./graph.h"

#include <tuple>
#include <vector>
#include <fstream>
#include <string>

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    void print_population();
    void evaluate_population();
    Solution& solve();

private:
    Model& model_ref;
    std::vector<Solution> population;
    Solution best_solution;
    Solution worst_solution;
    std::size_t avg_fitness;

    std::size_t fitness_evaluations {0};
    std::size_t generation_number {0};
    std::size_t population_size {500};
    std::size_t subgroup_size {7};
    double crossing_probability {0.8};
    double mutation_probability {0.1};

    std::string print_generation_info();
    bool check_reached_optimum();
    bool check_reached_ffe_limit();
    bool check_reached_gen_limit();

    void initialize_population(PopulationType population_type, std::size_t size);
    void random_initialization();
    void greedy_initialization();
    void mixed_initialization();

    std::vector<Solution> tournament_selection(std::size_t subgroup_size);
    std::vector<Solution> process_crossover(const std::vector<Solution>& parents);
    void normalize_parent_colours(const std::vector<Solution*>& parents);
    void evolve_population(
        std::vector<Solution>& parents, std::vector<Solution>& offsprings);
    void process_mutation();

    Solution double_point_crossover(
        const Graph& first_parent_graph,
        const Graph& second_parent_graph,
        const std::size_t first_crossing_point,
        const std::size_t second_crossing_point);
    Solution create_new_solution(Graph&& graph);
    double fitness_evaluation(Solution& solution);
};
