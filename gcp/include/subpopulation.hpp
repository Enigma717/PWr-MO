#pragma once

#include "structs/solution.hpp"

#include <vector>

class Model;

class Subpopulation {
public:
    using BuildingBlocks = std::vector<std::vector<std::size_t>>;

    Subpopulation() = delete;
    Subpopulation(std::size_t subpopulation_size, Model& model_ref);

    std::size_t subpopulation_size {2uz};
    std::size_t iterations_done {0uz};
    bool is_locked {false};

    std::vector<Solution> individuals;
    std::vector<Solution> improving_offsprings;
    Solution* best_solution;
    Solution* worst_solution;
    double avg_fitness {0.0};

    static std::size_t fitness_evaluations;
    static std::size_t get_ffe();

    void print_individuals() const;
    void print_info() const;
    void run_iteration();

private:
    Model& model_ref;

    Solution create_new_solution(Graph&& graph);
    double fitness_evaluation(Solution& solution);
    void subsitute_subpopulation_with_offsprings();
    void update_subpopulation_data();

    std::vector<Solution*> tournament_selection();
    void process_crossover(std::vector<Solution*> candidates);
    void process_partition_crossover(const std::vector<Solution*>& parents);
    BuildingBlocks normalize_parent_colours(
        Graph& first_parent_graph,
        Graph& second_parent_graph);
};
