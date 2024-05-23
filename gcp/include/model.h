#pragma once

#include "structs/model_params.h"
#include "loader.h"
#include "solver.h"

#include <memory>
#include <random>

class Graph;

class Model {
public:
    Model();

    void load_file(const std::string& file_path);
    void create_graph();
    void add_edge_to_graph(
        const std::size_t source_id,
        const std::size_t destination_id);

    std::string print_model_parms() const;

    bool check_colouring_corretness() const;
    std::size_t evaluate_fitness() const;

    void solve_random();

public:
    ModelParams model_params;
    Loader loader;
    std::unique_ptr<Graph> graph;
    Solver solver;

    std::random_device rd;
    std::mt19937 rng;

    friend class Loader;
};
