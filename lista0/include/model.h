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
    Model();

    void load_file(const std::string& file_path);
    void create_weight_matrix();
    void solve_knapsack();

    std::string print_model_parms() const;
    std::string print_nodes() const;
    std::string print_items() const;
    std::string print_weight_matrix() const;

    double objective_function(const std::vector<Node>& solution) const;
    double evaluate_member_fitness(const Member& member);
    double calculate_travel_time(const double distance, const int item_weight);

    std::vector<Node> k_random_solution(std::size_t k_factor);
    std::vector<Node> nearest_neighbour(std::size_t starting_point) const;
    std::vector<Node> extended_nearest_neighbour() const;

private:
    ModelParams model_params;
    Loader loader;

    std::random_device rd;
    std::mt19937 rng{rd()};
    std::vector<Node> nodes;
    std::vector<Item> items;
    std::vector<std::vector<double>> weights;
    // std::optional<double> optimum {std::nullopt};

    friend class Loader;
    friend class GeneticSolver;

public:
    std::vector<bool> knapsack_solution;
    int knapsack_value {0};
    int current_knapsack_weight {0};

    GeneticSolver genetic_solver;
};
