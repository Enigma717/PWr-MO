#pragma once

#include "structs/model_params.h"
#include "structs/solution.h"
#include "genetic_solver.h"
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
    void create_base_graph();
    void add_edge_to_base_graph(
        const std::size_t source_id,
        const std::size_t destination_id);

    std::size_t calculate_max_degree() const;
    std::string print_model_parms() const;

    bool check_colouring_corretness(const Graph& solution) const;
    bool check_no_colouring_collisions(const std::vector<std::size_t>& forbidden_colours, std::size_t vertex_colour) const;
    std::set<std::size_t> get_used_colours(const Graph& solution) const;
    std::vector<std::size_t> get_forbidden_colours(const Vertex& vertex) const;
    std::size_t find_available_colour(const std::vector<std::size_t>& forbidden_colours) const;
    Graph& fix_colouring(Graph& solution) const;
    std::size_t evaluate_fitness(Graph& solution) const;

    Graph solve_random();
    Graph solve_random(Graph graph);
    Graph solve_greedy();
    Graph solve_greedy(Graph graph);
    Solution& solve_genetic();

public:
    ModelParams model_params;
    Loader loader;
    std::unique_ptr<Graph> base_graph;
    Solver solver;
    GeneticSolver genetic_solver;

    std::random_device rd;
    std::mt19937 rng;
    std::set<std::size_t> final_colours;

    friend class Loader;
};
