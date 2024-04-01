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
    constexpr std::size_t generation_limit {500'000uz};
    constexpr std::size_t ffe_limit {5'000'000uz};
    constexpr std::uint64_t k_random_factor {100u};
    constexpr double random_initialization_probability {0.9};

    constexpr double ln_factor(const std::size_t iteration, const std::size_t dimension)
    {
        return (std::log(iteration) / std::log(dimension));
    }
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
    const Member new_best_member {*std::max_element(population.begin(), population.end())};
    const Member worst_member {*std::min_element(population.begin(), population.end())};
    const double avg_fitness {
        std::accumulate(population.begin(), population.end(), 0.0) / population_size};

    if (generation_number == 0)
        best_member = new_best_member;

    if (new_best_member > best_member) {
        best_member = new_best_member;

        std::cout << std::fixed << std::setprecision(2)
            << "|-> Generation number: " << generation_number
            << "\t||\tBest fitness: " << best_member.fitness
            << "\tKnapsack value: " << best_member.knapsack_value
            << "\tTraveling time: " << best_member.knapsack_value - best_member.fitness
            << "\t||\tWorst fitness: " << worst_member.fitness
            << "\t||\tAverage fitness: " << avg_fitness
            << "\t||\tFFE: " << fitness_evaluations << "\n";
    }
}

const Member& GeneticSolver::solve()
{
    initialize_population(PopulationType::mixed, population_size);
    evaluate_population();

    std::cout << "\n===========================\n\n";

    while (fitness_evaluations < ffe_limit && generation_number < generation_limit) {
        const std::vector<Member> parents {tournament_selection(subgroup_size)};
        const std::vector<Member> offsprings {process_crossover({parents[0], parents[1]})};

        evolve_population(parents, offsprings);
        process_mutation();
        evaluate_population();

        generation_number++;
    }

    std::cout << "\\-> Last generation number: " << generation_number << "\n";

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
                member.route.begin() + first_inverse_point,
                member.route.begin() + second_inverse_point);

            member = create_new_member(member.route);
        }
    }
}


void GeneticSolver::random_initialization()
{
    for (auto& member : population)
        member = create_new_member(model_ref.k_random_solution(k_random_factor));
}

void GeneticSolver::neighbour_initialization()
{
    for (std::size_t member_position {0uz}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};
        const std::uint16_t starting_node_index {
            static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

        member = create_new_member(model_ref.nearest_neighbour(starting_node_index));
    }
}

void GeneticSolver::mixed_initialization()
{
    for (std::size_t member_position {0uz}; member_position < population.size(); member_position++) {
        auto& member {population.at(member_position)};

        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        const double random_choice {distribution(model_ref.rng)};

        if (random_choice < random_initialization_probability) {
            member = create_new_member(model_ref.k_random_solution(k_random_factor));
        }
        else {
            const std::uint16_t starting_node_index {
                static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

            member = create_new_member(model_ref.nearest_neighbour(starting_node_index));
        }
    }
}

Member GeneticSolver::order_crossover(
    const Member& first_parent,
    const Member& second_parent,
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

    return create_new_member(offspring_solution);
}

Member GeneticSolver::create_new_member(const std::vector<Node>& route)
{
    Member member;
    member.route = route;
    member.penalized_items = penalize_item_values(route);
    member.packing_plan = model_ref.solve_knapsack_greedy(member);
    member.fitness = fitness_evaluation(member);

    return member;
}

std::vector<Item> GeneticSolver::penalize_item_values(const std::vector<Node>& route)
{
    std::vector<Item> penalized_items {model_ref.items};

    for (std::size_t i {0uz}; i < route.size(); i++) {
        const std::size_t current_node {route.at(i).index};

        if (current_node != 1) {
            Item& current_item {penalized_items.at(current_node - 2)};
            current_item.profit *= ln_factor(i + 2, route.size() + 1);
            current_item.ratio =
                static_cast<double>(current_item.profit) /
                static_cast<double>(current_item.weight);
        }
    }

    return penalized_items;
}


double GeneticSolver::fitness_evaluation(Member& member)
{
    fitness_evaluations++;
    return model_ref.evaluate_member_fitness(member);
}
