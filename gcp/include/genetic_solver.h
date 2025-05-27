#pragma once

#include "./enums/crossover_type.h"
#include "./enums/population_type.h"
#include "./structs/solution.h"
#include "./graph.h"

#include <vector>
#include <fstream>
#include <string>

using BuildingBlocks = std::vector<std::vector<std::size_t>>;

class Model;

class GeneticSolver {
public:
    GeneticSolver() = delete;
    GeneticSolver(Model& model_ref);

    double variance();
    Solution& solve(double& avg);

    std::size_t fitness_evaluations {0uz};
    std::size_t generation_number {0uz};
    std::size_t population_size {200uz};
    std::size_t tournament_size {4uz};
    double crossing_probability {0.8};
    double mutation_probability {0.2};
    CrossoverType crossover_type {CrossoverType::partition};

private:
    Model& model_ref;
    std::vector<Solution> population;
    Solution* best_solution;
    Solution* worst_solution;
    double avg_fitness;

    std::string print_generation_info();
    bool check_reached_optimum();
    bool check_reached_ffe_limit();
    bool check_reached_gen_limit();

    void initialize_population(PopulationType population_type, std::size_t size);
    void random_initialization();
    void greedy_initialization();
    void mixed_initialization();

    void print_population();
    void evaluate_population(std::ofstream& csv_file);

    std::vector<Solution> tournament_selection();
    std::vector<Solution> crossover_parents(std::vector<Solution>& parents);
    std::vector<Solution> process_crossover(
        Graph& first_parent_graph, Graph& second_parent_graph);
    void evolve_population(
        std::vector<Solution>& parents, std::vector<Solution>& offsprings);
    void process_mutation();

    std::vector<Solution> process_double_point_crossover(
        const Graph& first_parent_graph, const Graph& second_parent_graph);
    Solution double_point_crossover(
        const Graph& first_parent_graph,
        const Graph& second_parent_graph,
        const std::size_t first_crossing_point,
        const std::size_t second_crossing_point);

    std::vector<Solution> process_uniform_crossover(
        const Graph& first_parent_graph, const Graph& second_parent_graph);
    Solution uniform_crossover(
        const Graph& first_parent_graph, const Graph& second_parent_graph);

    std::vector<Solution> process_partition_crossover(
        Graph& first_parent_graph, Graph& second_parent_graph);
    BuildingBlocks normalize_parent_colours(
        Graph& first_parent_graph, Graph& second_parent_graph);

    Solution create_new_solution(Graph&& graph);
    double fitness_evaluation(Solution& solution);
};
