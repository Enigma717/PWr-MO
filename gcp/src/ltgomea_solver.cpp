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
#include <unordered_map>
#include <utility>

namespace
{
    constexpr std::size_t subpopulation_base_size {2uz};
    constexpr std::size_t ltgomea_base_iterations {4uz};
    constexpr std::size_t iterations_limit {100'000uz};
    constexpr std::size_t ffe_limit {1'00'000uz};
}

LTGOMEASolver::LTGOMEASolver(Model& model_ref) : model_ref{model_ref} {}

void LTGOMEASolver::solve()
{
    std::string dir_type {crossover_type == CrossoverType::optimal_mixing ? "ltgomea" : "ltgomeapx"};

    std::ofstream results_file;
    std::stringstream results_path;
    results_path << "./csv/results/" + dir_type + "/results_" << model_ref.model_params.instance_name << ".csv";
    results_file.open(results_path.str(), std::ios_base::app);

    if(!results_file.is_open()) {
        perror("Error opening results file");
        exit(EXIT_FAILURE);
    }

    std::ofstream plot_file;
    std::stringstream plot_data_path;
    plot_data_path << "./csv/results/" + dir_type + "/plot_" << model_ref.model_params.instance_name << ".csv";
    plot_file.open(plot_data_path.str(), std::ios_base::app);

    if(!plot_file.is_open()) {
        perror("Error opening plot file");
        exit(EXIT_FAILURE);
    }

    plot_file << "iteration; subpops; biggest_subpop; best; avg; dev\n";
    results_file << "subpops; biggest_subpop; best; avg; dev; iterations; ffe; time\n";

    total_iterations = 0uz;
    Subpopulation::fitness_evaluations = 0;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    while (!is_optimum_reached
        && total_iterations < iterations_limit
        && Subpopulation::get_ffe() < ffe_limit
    ) {
        create_new_subpopulation();
        generational_step(subpopulations_count - 1, plot_file);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    results_file << subpopulations_count << "; "
        << subpopulations.back().subpopulation_size << "; "
        << best_solution->fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "; "
        << total_iterations << "; "
        << Subpopulation::get_ffe() << "; "
        << elapsed_time << "\n";

    plot_file << total_iterations << "; "
        << subpopulations_count << "; "
        << subpopulations.back().subpopulation_size << "; "
        << best_solution->fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "\n";

    plot_file.close();
    results_file.close();

    std::cout << "\nFinal solution: [" << best_solution << "] | " << *best_solution;
    std::cout << "\nUsed crossover type: [" << (crossover_type == CrossoverType::optimal_mixing ? "optimal_mixing]\n" : "partition_crossover]\n");

    return;
}

void LTGOMEASolver::calculate_solver_info()
{
    avg_fitness = 0.0;
    best_solution = subpopulations.at(0).best_solution;

    for (const auto& subpopulation : subpopulations) {
        if (subpopulation.best_solution->fitness < best_solution->fitness)
            best_solution = subpopulation.best_solution;

        avg_fitness += subpopulation.avg_fitness;

        double sum {0.0};

        for (const auto& solution : subpopulation.individuals)
            sum += std::pow((solution.fitness - subpopulation.avg_fitness), 2);

        avg_deviation += sqrt(sum / subpopulation.subpopulation_size);
    }

    avg_fitness /= subpopulations_count;
    avg_deviation /= subpopulations_count;
}

void LTGOMEASolver::create_new_subpopulation()
{
    if (subpopulations_count == 0)
        subpopulations.emplace_back(
            subpopulation_base_size, crossover_type, model_ref);
    else
        subpopulations.emplace_back(
            subpopulation_base_size * subpopulations.back().subpopulation_size,
            crossover_type,
            model_ref);

    subpopulations_count++;
}

void LTGOMEASolver::generational_step(std::size_t current_subpopulation_index, std::ofstream& plot_file)
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

            calculate_solver_info();

            plot_file << total_iterations << "; "
                << subpopulations_count << "; "
                << subpopulations.back().subpopulation_size << "; "
                << best_solution->fitness << "; "
                << avg_fitness << "; "
                << avg_deviation << "\n";

            current_subpopulation.print_info();

            total_iterations++;
        }

        if (current_subpopulation_index > 0)
            generational_step(current_subpopulation_index - 1, plot_file);
    }
}

bool LTGOMEASolver::check_stop_condition(std::size_t current_subpopulation_index)
{
    const auto& current_subpopulation {subpopulations.at(current_subpopulation_index)};

    for (std::size_t larger_subpopulations_index = current_subpopulation_index + 1;
        larger_subpopulations_index < subpopulations_count;
        larger_subpopulations_index++
    ) {
        if (current_subpopulation.avg_fitness
            > subpopulations.at(larger_subpopulations_index).avg_fitness
        )
            return true;
    }

    for (std::size_t i {1uz}; i < current_subpopulation.individuals.size(); i++) {
        if (model_ref.are_two_solutions_same(
            current_subpopulation.individuals.at(0), current_subpopulation.individuals.at(i))
        )
            return true;
    }

    return false;
}
