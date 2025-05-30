#include "genetic_solver.hpp"
#include "model.hpp"
#include "utility_operators.hpp"

#include <iostream>

namespace
{
    constexpr std::size_t subpopulation_base_size {2uz};
    constexpr std::size_t ltgomea_base_iterations {4uz};
    constexpr std::size_t iterations_limit {1'000'000uz};
    constexpr std::size_t ffe_limit {500'000uz};
}

LTGOMEASolver::LTGOMEASolver(Model& model_ref) : model_ref{model_ref} {}

void LTGOMEASolver::solve()
{
    total_iterations = 0uz;

    while (!is_optimum_reached
        && total_iterations < iterations_limit
        && Subpopulation::get_ffe() < ffe_limit
    ) {
        create_new_subpopulation();

        generational_step(subpopulations_count - 1);
    }

    Solution* final_solution {subpopulations.at(0).best_solution};

    std::cout << "\n\nFinal results:\n";
    for (const auto& subpopulation : subpopulations) {
        if (subpopulation.best_solution->fitness < final_solution->fitness)
            final_solution = subpopulation.best_solution;

        subpopulation.print_info();
    }

    std::cout << "\nFinal solution: [" << final_solution << "] | " << *final_solution;

    return;
}

void LTGOMEASolver::create_new_subpopulation()
{
    if (subpopulations_count == 0)
        subpopulations.emplace_back(subpopulation_base_size, model_ref);
    else
        subpopulations.emplace_back(subpopulation_base_size * subpopulations.back().subpopulation_size, model_ref);

    subpopulations_count++;
}

void LTGOMEASolver::generational_step(std::size_t current_subpopulation_index)
{
    auto& current_subpopulation {subpopulations.at(current_subpopulation_index)};

    for (std::uint8_t i {0u}; i < ltgomea_base_iterations; i++) {
        if (!current_subpopulation.is_locked)
            current_subpopulation.is_locked = check_stop_condition(current_subpopulation_index);

        if (current_subpopulation.iterations_done > 0
            && current_subpopulation.best_solution->fitness == model_ref.model_params.optimum
        ) {
            is_optimum_reached = true;

            return;
        }

        if (total_iterations > iterations_limit)
            return;

        if (Subpopulation::get_ffe() > ffe_limit)
            return;

        if (!current_subpopulation.is_locked) {
            current_subpopulation.run_iteration();
            total_iterations++;
        }

        if (!current_subpopulation.is_locked) {
            current_subpopulation.print_info();
            std::cout << "Current FFE: " << Subpopulation::get_ffe() << "\n";
        }

        if (current_subpopulation_index > 0)
            generational_step(current_subpopulation_index - 1);
    }
}

bool LTGOMEASolver::check_stop_condition(std::size_t current_subpopulation_index)
{
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
