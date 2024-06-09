#pragma once

#include "structs/model_params.h"
#include "loader.h"
#include "solver.h"

#include <memory>
#include <random>
#include <set>

class Graph;

class Model {
public:
    Model();

    void load_file(const std::string& file_path);
    void create_graph();
    void add_edge_to_graph(
        const std::size_t source_id,
        const std::size_t destination_id);

    std::size_t calculate_max_degree() const;
    std::string print_model_parms() const;

    bool check_colouring_corretness(const Graph& solution) const;
    std::size_t evaluate_fitness(const Graph& solution);

    Graph solve_random(Graph graph);
    Graph solve_greedy(Graph graph);

public:
    ModelParams model_params;
    Loader loader;
    std::unique_ptr<Graph> graph;
    Solver solver;

    std::random_device rd;
    std::mt19937 rng;
    std::set<std::size_t> final_colours;

    friend class Loader;
};
