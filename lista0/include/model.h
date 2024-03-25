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

    std::string print_model_parms() const;
    std::string print_nodes() const;
    std::string print_weight_matrix() const;
    std::vector<Node> get_nodes() const;

    double objective_function(const std::vector<Node>& solution) const;
    double prd(double objective_sum) const;

    std::vector<Node> k_random_solution(std::uint64_t k);
    std::vector<Node> nearest_neighbour(std::uint16_t starting_point) const;
    std::vector<Node> extended_nearest_neighbour() const;

private:
    ModelParams model_params;
    Loader loader;

    std::random_device rd;
    std::mt19937 rng{rd()};
    std::vector<Node> nodes;
    std::vector<Item> items;
    std::vector<std::vector<double>> weights;

    friend class Loader;
    friend class GeneticSolver;

public:
    GeneticSolver genetic_solver;

};
