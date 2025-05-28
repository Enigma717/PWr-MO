#pragma once

#include "structs/solution.hpp"

#include <vector>

class Model;

class Subpopulation {
public:
    Subpopulation() = delete;
    Subpopulation(std::size_t subpopulation_size, Model& model_ref);

    std::size_t subpopulation_size {2uz};
    std::size_t fitness_evaluations {0uz};
    std::size_t iterations_done {0uz};
    double avg_fitness {0.0};
    bool locked {false};
    std::vector<Solution> individuals;

    void print_individuals() const;
    void print_info() const;
    void run_iteration();

private:
    Model& model_ref;
    Solution* best_solution;
    Solution* worst_solution;

    Solution create_new_solution(Graph&& graph);
    double fitness_evaluation(Solution& solution);
};
