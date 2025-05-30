#pragma once

#include "enums/crossover_type.hpp"
#include "enums/population_type.hpp"
#include "structs/solution.hpp"
#include "subpopulation.hpp"
#include "graph.hpp"

#include <vector>
#include <fstream>
#include <string>

class Model;

class LTGOMEASolver {
public:
    LTGOMEASolver() = delete;
    LTGOMEASolver(Model& model_ref);

    void solve();

    std::size_t total_iterations {0uz};
    std::size_t subpopulations_count {0uz};
    std::size_t subpopulations_limit {4uz};
    bool is_optimum_reached {false};

private:
    Model& model_ref;
    std::vector<Subpopulation> subpopulations;
    Solution* best_solution;
    Solution* worst_solution;
    double avg_fitness;

    void create_new_subpopulation();
    void generational_step(std::size_t largest_subpopulation_index);
    bool check_stop_condition(std::size_t current_subpopulation_index);
};
