#include "../include/genetic_solver.h"
#include "../include/model.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <random>
#include <ranges>
#include <algorithm>
#include <cmath>
#include <iomanip>

namespace
{
    constexpr std::size_t expected_winners {2uz};
    constexpr std::size_t dimension_divisor {4uz};
    constexpr std::size_t parents_pair_step {2uz};
    constexpr std::size_t generation_limit {50'000uz};
    constexpr std::size_t ffe_limit {500'000uz};
    constexpr std::uint64_t k_random_factor {100u};
    constexpr double random_initialization_probability {0.9};

    constexpr double variance(const std::vector<double>& values, const double mean)
    {
        double sum {0.0};

        for (const auto val : values) {
            sum += std::pow((val - mean), 2);
        }

        return (sum / 10.0);
    }
}

GeneticSolver::GeneticSolver(Model& model_ref) : model_ref{model_ref} {}

void GeneticSolver::print_population()
{
    for (std::size_t i {0uz}; i < population.size(); i++) {
        std:: cout << "Solution " << i << ": " << population.at(i) << "\n";
    }
}

void GeneticSolver::evaluate_population()
{
    const Solution new_best_solution {*std::max_element(population.begin(), population.end())};
    const Solution new_worst_solution {*std::min_element(population.begin(), population.end())};
    const double new_avg_fitness {
        std::accumulate(population.begin(), population.end(), 0.0) / population_size};

    if (generation_number == 0)
        best_solution = new_best_solution;

    if (new_best_solution > best_solution) {
        best_solution = new_best_solution;
        worst_solution = new_worst_solution;
        avg_fitness = new_avg_fitness;

        std::cout << std::fixed << std::setprecision(2)
            << "|-> Generation number: " << generation_number
            << "\t||\tBest fitness: " << best_solution.fitness
            << "\tKnapsack value: " << best_solution.knapsack_value
            << "\tTraveling time: " << best_solution.knapsack_value - best_solution.fitness
            << " \t||\tWorst fitness: " << worst_solution.fitness
            << "\t||\tAverage fitness: " << avg_fitness
            << "\t||\tFFE: " << fitness_evaluations << "\n";
    }
}

const Solution& GeneticSolver::solve()
{
    initialize_population(PopulationType::mixed, population_size);
    evaluate_population();

    std::cout << "\n===========================\n\n";

    while (fitness_evaluations < ffe_limit && generation_number < generation_limit) {
        const std::vector<Solution> parents {tournament_selection(subgroup_size)};
        const std::vector<Solution> offsprings {process_crossover(parents)};

        evolve_population(parents, offsprings);
        process_mutation();
        evaluate_population();

        generation_number++;
    }

    return best_solution;
}


void GeneticSolver::initialize_population(
    PopulationType population_type, std::size_t population_size)
{
    population.resize(population_size);

    switch (population_type) {
        case PopulationType::random: random_initialization(); break;
        case PopulationType::neighbour: neighbour_initialization(); break;
        case PopulationType::mixed: mixed_initialization(); break;
    }
}

void GeneticSolver::random_initialization()
{
    for (auto& solution : population)
        solution = create_new_solution(model_ref.k_random_solution(k_random_factor).route);
}

void GeneticSolver::neighbour_initialization()
{
    for (std::size_t solution_position {0uz}; solution_position < population.size(); solution_position++) {
        auto& solution {population.at(solution_position)};
        const std::uint16_t starting_node_index {
            static_cast<std::uint16_t>((solution_position % model_ref.model_params.dimension) + 1)};

        solution = create_new_solution(model_ref.nearest_neighbour(starting_node_index));
    }
}

void GeneticSolver::mixed_initialization()
{
    for (std::size_t solution_position {0uz}; solution_position < population.size(); solution_position++) {
        auto& solution {population.at(solution_position)};

        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        const double random_choice {distribution(model_ref.rng)};

        if (random_choice < random_initialization_probability) {
            solution = create_new_solution(model_ref.k_random_solution(k_random_factor).route);
        }
        else {
            const std::uint16_t starting_node_index {
                static_cast<std::uint16_t>((solution_position % model_ref.model_params.dimension) + 1)};

            solution = create_new_solution(model_ref.nearest_neighbour(starting_node_index));
        }
    }
}

std::vector<Solution> GeneticSolver::tournament_selection(std::size_t subgroup_size)
{
    const std::size_t subgroups_count {population.size() / subgroup_size};
    const int winners_count {static_cast<int>(subgroup_size - expected_winners)};
    const auto drop_point {std::max(0, winners_count)};
    std::vector<Solution> final_winners;

    final_winners.reserve(subgroups_count * expected_winners);

    for (const auto& subgroup : population | std::views::chunk(subgroup_size)) {
        std::sort(subgroup.begin(), subgroup.end());

        const auto& winners {std::views::drop(subgroup, drop_point)};

        for (const Solution& winner : winners)
            final_winners.push_back(winner);
    }

    return final_winners;
}

std::vector<Solution> GeneticSolver::process_crossover(const std::vector<Solution>& parents)
{
    const std::size_t dimension {model_ref.model_params.dimension};
    const std::size_t dimension_quarter {dimension / dimension_divisor};

    std::uniform_int_distribution<std::size_t> int_distribution(
        dimension_quarter, dimension - dimension_quarter);
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);

    std::vector<Solution> offsprings;

    offsprings.reserve(parents.size());

    for (std::size_t i {0uz}; i < parents.size(); i += parents_pair_step) {
        const double probability {real_distribution(model_ref.rng)};
        const Solution& first_parent {parents.at(i)};
        const Solution& second_parent {parents.at(i + 1)};

        if (probability < crossing_probability) {
            std::size_t first_crossing_point {int_distribution(model_ref.rng)};
            std::size_t second_crossing_point {int_distribution(model_ref.rng)};

            if (first_crossing_point > second_crossing_point)
                std::swap(first_crossing_point, second_crossing_point);

            const Solution first_offspring {order_crossover(
                first_parent, second_parent, first_crossing_point, second_crossing_point)};
            const Solution second_offspring {order_crossover(
                second_parent, first_parent, first_crossing_point, second_crossing_point)};

            offsprings.push_back(first_offspring);
            offsprings.push_back(second_offspring);
        }
        else {
            offsprings.push_back(first_parent);
            offsprings.push_back(second_parent);
        }
    }

    return offsprings;
}


void GeneticSolver::evolve_population(
    const std::vector<Solution>& parents, const std::vector<Solution>& offsprings)
{
    std::size_t current_position {0uz};

    while (current_position < parents.size()) {
        population.at(current_position) = parents.at(current_position);
        current_position++;
    }

    const std::size_t remaining_places {model_ref.model_params.dimension - current_position};
    std::size_t offspring_position {0uz};

    while (offspring_position < offsprings.size()) {
        if (offspring_position >= remaining_places)
            return;

        population.at(current_position) = offsprings.at(offspring_position);
        current_position++;
        offspring_position++;
    }
}

void GeneticSolver::process_mutation()
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (auto& solution : population) {
        const double probability {distribution(model_ref.rng)};

        if (probability < mutation_probability) {
            solution.route = model_ref.process_invert_mutation(solution.route);
            solution = create_new_solution(solution.route);
        }
    }
}

Solution GeneticSolver::order_crossover(
    const Solution& first_parent,
    const Solution& second_parent,
    const std::size_t first_crossing_point,
    const std::size_t second_crossing_point)
{
    const std::size_t dimension {first_parent.route.size()};
    std::vector<bool> used_nodes(dimension);
    std::vector<Node> offspring_solution(dimension);

    std::size_t current_position {first_crossing_point};

    while (current_position <= second_crossing_point) {
        offspring_solution.at(current_position) = first_parent.route.at(current_position);
        used_nodes.at(first_parent.route.at(current_position).index - 1) = true;

        current_position++;
    }

    std::size_t first_legal_position {0uz};
    current_position = 0uz;

    while (current_position < dimension) {
        if (current_position >= first_crossing_point && current_position <= second_crossing_point) {
            current_position++;
            continue;
        }

        while (used_nodes.at(second_parent.route.at(first_legal_position).index - 1))
            first_legal_position++;

        offspring_solution.at(current_position) = second_parent.route.at(first_legal_position);
        used_nodes.at(second_parent.route.at(first_legal_position).index - 1) = true;

        current_position++;
    }

    return create_new_solution(offspring_solution);
}

Solution GeneticSolver::create_new_solution(const std::vector<Node>& route)
{
    Solution solution;
    solution.route = route;
    solution.penalized_items = model_ref.penalize_item_values(route);
    solution.packing_plan = model_ref.solve_knapsack_greedy(solution);
    solution.fitness = fitness_evaluation(solution);

    return solution;
}

double GeneticSolver::fitness_evaluation(Solution& solution)
{
    fitness_evaluations++;
    return model_ref.evaluate_solution_fitness(solution);
}
