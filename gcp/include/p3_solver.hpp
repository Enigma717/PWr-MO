#pragma once

#include "structs/solution.hpp"
#include "utility_operators.hpp"

#include <vector>
#include <unordered_map>

class Model;

class P3Solver {
public:
    P3Solver() = delete;
    P3Solver(Model& model_ref);

    void print_info() const;
    void solve();

    std::size_t fitness_evaluations {0uz};
    std::size_t total_iterations {0uz};
    std::size_t pyramid_levels {0uz};
    bool is_optimum_reached {false};
    std::uint8_t crossover_type {0u};

private:
    Model& model_ref;
    std::vector<std::vector<Solution>> pyramid;
    std::vector<LinkageTreeBuilder> linkage_trees;
    std::unordered_map<std::size_t, std::size_t> known_solutions;

    Solution* best_solution;

    void create_new_level();
    void next_iteration();

    Solution create_new_solution(Graph&& graph);
    void add_solution_to_level(const Solution& solution, const std::size_t level);
    double fitness_evaluation(Solution& solution);
    void normalize_colours(Graph& first_parent_graph, Graph& second_parent_graph);

    void process_optimal_mixing(Solution& solution, std::size_t current_level);
};
