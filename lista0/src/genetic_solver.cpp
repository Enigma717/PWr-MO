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
        double random_choice {distribution(model_ref.rng)};

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
