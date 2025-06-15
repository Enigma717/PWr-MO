#pragma once

#include "enums/crossover_type.hpp"
#include "structs/solution.hpp"
#include "linkage_tree_builder.hpp"

#include <vector>
#include <map>

class Model;

class Subpopulation {
public:
    using BuildingBlocks = std::vector<std::vector<std::size_t>>;
    using ColoursMap = std::map<std::size_t, std::size_t>;
    using ColoursPair = std::pair<std::size_t, std::size_t>;

    Subpopulation() = delete;
    Subpopulation(
        std::size_t subpopulation_size, CrossoverType crossover_type, Model& model_ref);

    LinkageTreeBuilder lt_builder;
    CrossoverType crossover_type {CrossoverType::optimal_mixing};

    std::size_t subpopulation_size {2uz};
    std::size_t iterations_done {0uz};
    bool is_locked {false};

    std::vector<Solution> individuals;
    std::vector<Solution> offsprings;

    Solution* best_solution;
    Solution* worst_solution;
    double avg_fitness {0.0};
    double deviation {0.0};

    static std::size_t fitness_evaluations;
    static std::size_t get_ffe();

    void print_individuals() const;
    void print_info() const;
    void run_iteration();

    Solution create_new_solution(Graph&& graph);
    double fitness_evaluation(Solution& solution);
    void subsitute_subpopulation_with_offsprings();
    double calculate_deviation();
    void update_subpopulation_data();
    void normalize_colours(Graph& first_parent_graph, Graph& second_parent_graph);

    std::vector<Solution*> tournament_selection();
    void process_partition_crossover(const std::vector<Solution*>& candidates);
    void partition_crossover(const std::vector<Solution*>& parents);
    BuildingBlocks obtain_building_blocks(
        const Graph& first_parent_graph,
        const Graph& second_parent_graph);

    void process_optimal_mixing();

private:
    Model& model_ref;
};
