#pragma once

#include "./loader.h"
#include "./genetic_solver.h"
#include "./structs/model_params.h"
#include "./structs/node.h"
#include "./structs/item.h"

#include <vector>
#include <random>

class Model {
public:
    Model();

    void load_file(const std::string& file_path);
    void create_weight_matrix();

    std::string print_model_parms();
    std::string print_nodes();
    std::string print_nodes(const std::vector<Node>& solution);
    std::string print_weight_matrix();

    double objective_function();
    double objective_function(const std::vector<Node>& solution);
    double prd(double objective_sum);

    std::vector<Node> k_random_solution(std::uint64_t k);
    std::vector<Node> nearest_neighbour(std::uint16_t starting_point);
    std::vector<Node> extended_nearest_neighbour();

private:
    ModelParams model_params;
    Loader loader;
    GeneticSolver genetic_solver;

    std::random_device rd;
    std::mt19937 rng{rd()};
    std::vector<Node> nodes;
    std::vector<Item> items;
    std::vector<std::vector<double>> weights;

    friend class Loader;
    friend class GeneticSolver;
};
