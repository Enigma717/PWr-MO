#include "genetic_solver.hpp"
#include "model.hpp"
#include "utility_operators.hpp"

#include <iostream>
#include <chrono>
#include <random>
#include <ranges>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <utility>

namespace
{
    constexpr std::size_t subpopulation_base_size {2uz};
    constexpr std::size_t ltgomea_base_iterations {4uz};
    constexpr std::size_t parents_pair_step {2uz};
    constexpr std::size_t print_info_threshold {100uz};
    constexpr std::size_t generation_limit {100'000uz};
    constexpr std::size_t ffe_limit {100'000uz};
}

LTGOMEASolver::LTGOMEASolver(Model& model_ref) : model_ref{model_ref} {}

void LTGOMEASolver::solve()
{
    // while (fitness_evaluations < 3) {
        // create_new_subpopulation();
    // }

    for (int i = 0; i < 3; i++) {
        create_new_subpopulation();

        // for (const auto& subpopulation : subpopulations) {
            // subpopulations.back().print_individuals();
        // }

        generational_step(subpopulations_count - 1);
    }

    std::cout << "\n\nFinal results:\n";
    for (const auto& subpopulation : subpopulations) {
        // std::cout << "\n";
        // subpopulation.print_individuals();
        subpopulation.print_info();
    }

    return;
}

void LTGOMEASolver::create_new_subpopulation()
{
    // std::cout << "\n\nSubpopulation count: " << subpopulations_count << "\n";

    if (subpopulations_count == 0) {
        subpopulations.emplace_back(subpopulation_base_size, model_ref);
    }
    else {
        const auto last_gomea_size {subpopulations.back().subpopulation_size};
        // std::cout << "\nLast GOMEA size: " << last_gomea_size;
        // std::cout << "\nNew GOMEA size: " << (subpopulation_base_size * last_gomea_size) << "\n\n";

        subpopulations.emplace_back(subpopulation_base_size * last_gomea_size, model_ref);
    }

    subpopulations_count++;
}

void LTGOMEASolver::generational_step(std::size_t current_subpopulation_index)
{
    auto& current_subpopulation {subpopulations.at(current_subpopulation_index)};

    for (std::uint8_t i {0u}; i < ltgomea_base_iterations; i++) {
        if (!current_subpopulation.locked)
            current_subpopulation.locked = check_stop_condition(current_subpopulation_index);

        if (!current_subpopulation.locked)
            current_subpopulation.run_iteration();

        if (!current_subpopulation.locked)
            current_subpopulation.print_info();

        if (current_subpopulation_index > 0)
            generational_step(current_subpopulation_index - 1);
    }
}

bool LTGOMEASolver::check_stop_condition(std::size_t current_subpopulation_index)
{
    fitness_evaluations++;
    for (std::size_t larger_subpopulations_index = current_subpopulation_index + 1;
        larger_subpopulations_index < subpopulations_count;
        larger_subpopulations_index++
    ) {
        if (subpopulations.at(current_subpopulation_index).avg_fitness
            > subpopulations.at(larger_subpopulations_index).avg_fitness
        )
            return true;
    }

    return false;
}
