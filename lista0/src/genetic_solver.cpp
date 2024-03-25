#include "../include/genetic_solver.h"
#include "../include/model.h"
#include "../include/stream_operators.h"

#include <iostream>
#include <random>
#include <ranges>
#include <algorithm>

namespace
{
    const std::size_t expected_winners {2uz};
}

GeneticSolver::GeneticSolver(Model& model_ref) : model_ref{model_ref} {}

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

void GeneticSolver::print_population()
{
    for (std::size_t i {0uz}; i < population.size(); i++) {
        std:: cout << "Member " << i << ": " << population.at(i) << "\n";
    }
}

void GeneticSolver::random_initialization()
{
    for (auto& member : population) {
        member.solution = model_ref.k_random_solution(100);
        member.fitness = model_ref.objective_function(member.solution);
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

        const auto& winners {std::views::drop(std::views::reverse(subgroup), drop_point)};

        for (const Member& winner : winners)
            final_winners.push_back(winner);
    }

    for (const Member& winner : final_winners) {
        std::cout << "-> Final winner: " << winner << "\n";
    }
    std::cout << "\nSIZE: " << final_winners.size() << "\n";

    return final_winners;
}

OffspringsPair GeneticSolver::ordered_crossover(const Member& first_parent, const Member& second_parent)
// OffspringsPair GeneticSolver::ordered_crossover()
{
    const std::size_t dimension {model_ref.model_params.dimension};
    const std::size_t dimension_quarter {dimension / 4};

    std::vector<bool> used_nodes(dimension);
    std::vector<Node> first_offspring(dimension);
    std::vector<Node> second_offspring(dimension);

    std::uniform_int_distribution<std::size_t> distribution(
        dimension_quarter, dimension - dimension_quarter);
    const std::size_t crossing_point {distribution(model_ref.rng)};

    std::size_t i {0};

    for (;i < crossing_point; i++) {
        first_offspring[i] = first_parent.solution[i];
        second_offspring[i + crossing_point] = first_parent.solution[i];
        used_nodes[first_parent.solution[i].index - 1] = true;
    }

    std::size_t j {i};

    while (i < dimension) {
        while (used_nodes[second_parent.solution[j].index - 1] != 0) {
            if (j == dimension - 1)
                j = -1;

            j++;
        }


        second_offspring[i - crossing_point] = second_parent.solution[i];
        first_offspring[i] = second_parent.solution[i];
        used_nodes[second_parent.solution[j].index - 1] = true;

        if (j == dimension - 1)
            j = -1;

        j++;
        i++;
    }

    std::cout << "\n> Crossing point: " << crossing_point << "\n";
    std::cout << "> First parent: " << first_parent << "\n";
    std::cout << "> Second parent: " << second_parent << "\n";
    std::cout << "==========================\n";
    std::cout << "> First offspring: " << first_offspring << "\n";
    std::cout << "> Second offspring: " << second_offspring << "\n";

    return {first_offspring, second_offspring};
}

void GeneticSolver::neighbour_initialization()
{
    for (size_t member_position {0}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};
        const std::uint16_t starting_node_index {
            static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

        member.solution = model_ref.nearest_neighbour(starting_node_index);
        member.fitness = model_ref.objective_function(member.solution);
    }
}

void GeneticSolver::mixed_initialization()
{
    for (size_t member_position {0}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};

        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        const double random_choice {distribution(model_ref.rng)};

        if (random_choice < 0.9) {
            member.solution = model_ref.k_random_solution(100);
            member.fitness = model_ref.objective_function(member.solution);
        }
        else {
            const std::uint16_t starting_node_index {
                static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

            member.solution = model_ref.nearest_neighbour(starting_node_index);
            member.fitness = model_ref.objective_function(member.solution);
        }
    }
}
