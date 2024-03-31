#include "../include/genetic_solver.h"
#include "../include/model.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <random>
#include <ranges>
#include <algorithm>

namespace
{
    constexpr std::size_t expected_winners {2uz};
    constexpr std::size_t dimension_divisor {4uz};
    constexpr std::size_t parents_pair_step {2uz};
    constexpr std::size_t generation_limit {1'000'000uz};
    constexpr std::size_t ffe_limit {10'000'000uz};
    constexpr std::uint64_t k_random_factor {100u};
    constexpr double random_initialization_probability {0.9};
}

GeneticSolver::GeneticSolver(Model& model_ref) : model_ref{model_ref} {}

void GeneticSolver::print_population()
{
    for (std::size_t i {0uz}; i < population.size(); i++) {
        std:: cout << "Member " << i << ": " << population.at(i) << "\n";
    }
}

void GeneticSolver::evaluate_population()
{
    const auto new_best_member {*std::max_element(population.begin(), population.end())};

    if (generation_number == 0)
        best_member = new_best_member;

    if (new_best_member > best_member) {
        best_member = new_best_member;

        std::cout << "|-> Generation number: " << generation_number
            << "\tFitness: " << best_member.fitness
            << "\tFFE: " << fitness_evaluations << "\n";
    }
}

const Member& GeneticSolver::solve()
{
    initialize_population(PopulationType::mixed, population_size);
    evaluate_population();

    std::cout << "\n===========================\n\n";

    while (fitness_evaluations < ffe_limit) {
        const std::vector<Member> parents {tournament_selection(subgroup_size)};
        const std::vector<Member> offsprings {process_crossover({parents[0], parents[1]})};

        evolve_population(parents, offsprings);
        process_mutation();
        evaluate_population();

        generation_number++;
    }

    return best_member;
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

std::vector<Member> GeneticSolver::tournament_selection(std::size_t subgroup_size)
{
    const std::size_t subgroups_count {population.size() / subgroup_size};
    const int winners_count {static_cast<int>(subgroup_size - expected_winners)};
    const auto drop_point {std::max(0, winners_count)};
    std::vector<Member> final_winners;

    final_winners.reserve(subgroups_count * expected_winners);

    for (const auto& subgroup : population | std::views::chunk(subgroup_size)) {
        std::sort(subgroup.begin(), subgroup.end());

        const auto& winners {std::views::drop(subgroup, drop_point)};

        for (const Member& winner : winners)
            final_winners.push_back(winner);
    }

    return final_winners;
}

std::vector<Member> GeneticSolver::process_crossover(const std::vector<Member>& parents)
{
    const std::size_t dimension {model_ref.model_params.dimension};
    const std::size_t dimension_quarter {dimension / dimension_divisor};

    std::uniform_int_distribution<std::size_t> int_distribution(
        dimension_quarter, dimension - dimension_quarter);
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);

    std::vector<Member> offsprings;

    offsprings.reserve(parents.size());

    for (std::size_t i {0uz}; i < parents.size(); i += parents_pair_step) {
        const double probability {real_distribution(model_ref.rng)};
        const Member& first_parent {parents.at(i)};
        const Member& second_parent {parents.at(i + 1)};

        if (probability < crossing_probability) {
            std::size_t first_crossing_point {int_distribution(model_ref.rng)};
            std::size_t second_crossing_point {int_distribution(model_ref.rng)};

            if (first_crossing_point > second_crossing_point)
                std::swap(first_crossing_point, second_crossing_point);

            const Member first_offspring {order_crossover(
                first_parent, second_parent, first_crossing_point, second_crossing_point)};
            const Member second_offspring {order_crossover(
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
    const std::vector<Member>& parents, const std::vector<Member>& offsprings)
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
    std::uniform_int_distribution<std::size_t> int_distribution(
        0, model_ref.model_params.dimension);

    for (auto& member : population) {
        const double probability {distribution(model_ref.rng)};

        if (probability < mutation_probability) {
            std::size_t first_inverse_point {int_distribution(model_ref.rng)};
            std::size_t second_inverse_point {int_distribution(model_ref.rng)};

            if (first_inverse_point > second_inverse_point)
                std::swap(first_inverse_point, second_inverse_point);

            std::reverse(
                member.solution.begin() + first_inverse_point,
                member.solution.begin() + second_inverse_point);

            fitness_evaluation(member);
        }
    }
}


Member GeneticSolver::order_crossover(
    const Member& first_parent,
    const Member& second_parent,
    const std::size_t first_crossing_point,
    const std::size_t second_crossing_point)
{
    const std::size_t dimension {first_parent.solution.size()};
    std::vector<bool> used_nodes(dimension);
    std::vector<Node> offspring_solution(dimension);

    std::size_t current_position {first_crossing_point};

    while (current_position <= second_crossing_point) {
        offspring_solution.at(current_position) = first_parent.solution.at(current_position);
        used_nodes.at(first_parent.solution.at(current_position).index - 1) = true;

        current_position++;
    }

    std::size_t first_legal_position {0uz};
    current_position = 0uz;

    while (current_position < dimension) {
        if (current_position >= first_crossing_point && current_position <= second_crossing_point) {
            current_position++;
            continue;
        }

        while (used_nodes.at(second_parent.solution.at(first_legal_position).index - 1))
            first_legal_position++;

        offspring_solution.at(current_position) = second_parent.solution.at(first_legal_position);
        used_nodes.at(second_parent.solution.at(first_legal_position).index - 1) = true;

        current_position++;
    }

    Member offspring {offspring_solution, 0.0};
    fitness_evaluation(offspring);

    return offspring;
}

void GeneticSolver::random_initialization()
{
    for (auto& member : population) {
        member.solution = model_ref.k_random_solution(k_random_factor);
        member.fitness = model_ref.evaluate_member_fitness(member);
    }
}

void GeneticSolver::neighbour_initialization()
{
    for (std::size_t member_position {0uz}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};
        const std::uint16_t starting_node_index {
            static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

        member.solution = model_ref.nearest_neighbour(starting_node_index);
        member.fitness = model_ref.evaluate_member_fitness(member);
    }
}

void GeneticSolver::mixed_initialization()
{
    for (std::size_t member_position {0uz}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};

        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        const double random_choice {distribution(model_ref.rng)};

        if (random_choice < random_initialization_probability) {
            member.solution = model_ref.k_random_solution(k_random_factor);
            member.fitness = model_ref.evaluate_member_fitness(member);
        }
        else {
            const std::uint16_t starting_node_index {
                static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

            member.solution = model_ref.nearest_neighbour(starting_node_index);
            member.fitness = model_ref.evaluate_member_fitness(member);
        }
    }
}

void GeneticSolver::fitness_evaluation(Member& member)
{
    member.fitness = model_ref.evaluate_member_fitness(member);
    fitness_evaluations++;
}
