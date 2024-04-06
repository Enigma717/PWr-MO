#pragma once

#include "./loader.h"
#include "./genetic_solver.h"
#include "./structs/model_params.h"
#include "./structs/node.h"
#include "./structs/item.h"

#include <vector>
#include <optional>
#include <random>

class Model {
public:
    using PackingPlan = std::vector<bool>;

    Model();

    void load_file(const std::string& file_path);
    void create_weight_matrix();

    std::string print_model_parms() const;
    std::string print_nodes() const;
    std::string print_items() const;
    std::string print_weight_matrix() const;

    double objective_function(const std::vector<Node>& route) const;
    double evaluate_solution_fitness(Solution& solution);

    Solution k_random_solution(std::size_t k_factor);
    std::vector<Node> nearest_neighbour(std::size_t starting_point);
    Solution extended_nearest_neighbour();
    Solution simulated_annealing();
    Solution genetic_algorithm();

private:
    ModelParams model_params;
    Loader loader;
    GeneticSolver genetic_solver;

    std::random_device rd;
    std::mt19937 rng{rd()};
    std::vector<Node> nodes;
    std::vector<Item> items;
    std::vector<std::vector<double>> weights;

    PackingPlan solve_knapsack_greedy(Solution& solution);
    PackingPlan solve_knapsack_dp(Solution& solution);

    double calculate_travel_time(
        const double distance, const int item_weight, std::size_t& knapsack_weight);

    std::vector<Node> process_invert_mutation(const std::vector<Node>& route);
    std::vector<Item> penalize_item_values(const std::vector<Node>& route);

    friend class Loader;
    friend class GeneticSolver;
};
