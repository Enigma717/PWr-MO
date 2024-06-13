#include "../include/genetic_solver.h"
#include "../include/model.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <random>
#include <ranges>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <tuple>

namespace
{
    constexpr std::size_t expected_winners {2uz};
    constexpr std::size_t dimension_divisor {4uz};
    constexpr std::size_t parents_pair_step {2uz};
    constexpr std::size_t print_info_threshold {100uz};
    constexpr std::size_t generation_limit {50'000uz};
    constexpr std::size_t ffe_limit {500'000uz};
    constexpr double random_initialization_probability {0.9};
}

GeneticSolver::GeneticSolver(Model& model_ref) : model_ref{model_ref} {}

void GeneticSolver::print_population()
{
    for (std::size_t i {0uz}; i < population.size(); i++) {
        std::cout << "Solution " << i << " (" << &population.at(i) << "): " << population.at(i) << "\n";
    }
}

std::string GeneticSolver::print_generation_info()
{
    std::stringstream log;

    log << std::fixed << std::setprecision(2)
        << "|-> Generation number: " << generation_number
        << "\t||\tBest fitness: " << best_solution.fitness
        << " \t||\tWorst fitness: " << worst_solution.fitness
        << "\t||\tAverage fitness: " << avg_fitness
        << "\t||\tFFE: " << fitness_evaluations << "\n";

    return log.str();
}

bool GeneticSolver::check_reached_optimum()
{
    return best_solution.fitness == model_ref.model_params.optimum;
}

bool GeneticSolver::check_reached_ffe_limit()
{
    return fitness_evaluations >= ffe_limit;
}

bool GeneticSolver::check_reached_gen_limit()
{
    return generation_number >= generation_limit;
}

void GeneticSolver::evaluate_population()
{
    auto new_best_solution {*std::min_element(population.begin(), population.end())};
    auto new_worst_solution {*std::max_element(population.begin(), population.end())};
    std::size_t new_avg_fitness {
        std::accumulate(population.begin(), population.end(), 0.0) / population_size};

    if (generation_number == 0)
        best_solution = new_best_solution;

    if (new_best_solution < best_solution) {
        std::swap(best_solution, new_best_solution);
        std::swap(worst_solution, new_worst_solution);
        std::swap(avg_fitness, new_avg_fitness);
        std::cout << print_generation_info();
    }

    if (generation_number % print_info_threshold == 0 && generation_number != 0)
        std::cout << print_generation_info();
}

Solution& GeneticSolver::solve()
{
    initialize_population(PopulationType::mixed, population_size);
    evaluate_population();
    // print_population();

    std::cout << "\n\n===========================\n\n";

    while (!check_reached_ffe_limit() && !check_reached_gen_limit() && !check_reached_optimum()) {
        std::vector<Solution> parents {tournament_selection(subgroup_size)};
        std::vector<Solution> offsprings {process_crossover(parents)};

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
        case PopulationType::greedy: greedy_initialization(); break;
        case PopulationType::mixed: mixed_initialization(); break;
    }
}

void GeneticSolver::random_initialization()
{
    for (auto& solution : population)
        solution = create_new_solution(model_ref.solve_random());
}

void GeneticSolver::greedy_initialization()
{
    for (auto& solution : population)
        solution = create_new_solution(model_ref.solve_greedy());
}

void GeneticSolver::mixed_initialization()
{
    for (std::size_t solution_position {0uz}; solution_position < population.size(); solution_position++) {
        auto& solution {population.at(solution_position)};

        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        const double random_choice {distribution(model_ref.rng)};

        if (random_choice < random_initialization_probability)
            solution = create_new_solution(model_ref.solve_random());
        else
            solution = create_new_solution(model_ref.solve_greedy());
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
        std::sort(subgroup.begin(), subgroup.end(), std::greater{});

        const auto& winners {std::views::drop(subgroup, drop_point)};

        for (const Solution& winner : winners)
            final_winners.push_back(winner);
    }

    return final_winners;
}

std::vector<Solution> GeneticSolver::process_crossover(const std::vector<Solution>& parents)
{
    const std::size_t dimension {model_ref.model_params.vertices};
    const std::size_t dimension_quarter {dimension / dimension_divisor};

    std::uniform_int_distribution<std::size_t> int_distribution(
        dimension_quarter, dimension - dimension_quarter);
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);

    std::vector<Solution> offsprings;
    offsprings.reserve(parents.size());

    for (std::size_t i {0uz}; i < parents.size(); i += parents_pair_step) {
        const double probability {real_distribution(model_ref.rng)};
        const Solution first_parent {parents.at(i)};
        const Solution second_parent {parents.at(i + 1)};
        const Graph& first_parent_graph {first_parent.graph};
        const Graph& second_parent_graph {second_parent.graph};

        if (probability < crossing_probability) {
            std::size_t first_crossing_point {int_distribution(model_ref.rng)};
            std::size_t second_crossing_point {int_distribution(model_ref.rng)};

            if (first_crossing_point > second_crossing_point)
                std::swap(first_crossing_point, second_crossing_point);

            const Solution first_offspring {double_point_crossover(
                first_parent_graph, second_parent_graph, first_crossing_point, second_crossing_point)};
            const Solution second_offspring {double_point_crossover(
                second_parent_graph, first_parent_graph, first_crossing_point, second_crossing_point)};

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


Solution GeneticSolver::double_point_crossover(
    const Graph& first_parent_graph,
    const Graph& second_parent_graph,
    const std::size_t first_crossing_point,
    const std::size_t second_crossing_point)
{
    const std::size_t dimension {first_parent_graph.vertices.size()};
    Graph offspring_graph {first_parent_graph};
    std::size_t current_position {0uz};

    while (current_position < dimension) {
        std::size_t colour;

        if (current_position >= first_crossing_point && current_position <= second_crossing_point)
            colour = *second_parent_graph.colours.at(current_position);
        else
            colour = *first_parent_graph.colours.at(current_position);

        offspring_graph.vertices.at(current_position).update_colour(colour);

        current_position++;
    }

    return create_new_solution(std::move(offspring_graph));
}

// void GeneticSolver::normalize_parent_colours(const std::vector<Solution*>& parents)
// {
//     const Solution* const first_parent {parents.at(0)};
//     const Solution* const second_parent {parents.at(1)};
//     const Graph& first_parent_graph {first_parent->graph};
//     const Graph& second_parent_graph {second_parent->graph};

//     std::cout << "\n==============[NORMALIZATION]=============\n\n";
//     std::cout << "\nFirst parent (" << first_parent << "): " << *first_parent;
//     std::cout << "\nSecond parent (" << second_parent << "): " << *second_parent;

//     using ColoursMap = std::map<std::size_t, std::size_t>;
//     using ColoursPair = std::tuple<std::size_t, std::size_t>;

//     ColoursMap first_parent_colour_occurrences;
//     ColoursMap second_parent_colour_occurrences;
//     std::map<ColoursPair, std::size_t> common_colours;

//     for (std::size_t i {0uz}; i < first_parent_graph.vertices.size(); i++)
//     {
//         std::size_t colour {*first_parent_graph.colours.at(i)};

//         auto it(first_parent_colour_occurrences.find(colour));
//         if (it != first_parent_colour_occurrences.end())
//            it->second++;
//         else
//            first_parent_colour_occurrences[colour] = 1;

//     }

//     for (std::size_t i {0uz}; i < second_parent_graph.vertices.size(); i++)
//     {
//         std::size_t colour {*second_parent_graph.colours.at(i)};

//         auto it(second_parent_colour_occurrences.find(colour));
//         if (it != second_parent_colour_occurrences.end())
//            it->second++;
//          else
//            second_parent_colour_occurrences[colour] = 1;

//     }

//     for (std::size_t i {0uz}; i < first_parent_graph.vertices.size(); i++)
//     {
//         std::size_t first_parent_colour {*first_parent_graph.colours.at(i)};
//         std::size_t second_parent_colour {*second_parent_graph.colours.at(i)};
//         ColoursPair colours_pair {std::make_tuple(first_parent_colour, second_parent_colour)};

//         auto it(common_colours.find(colours_pair));

//         if (it != common_colours.end())
//            it->second++;
//         else
//            common_colours[colours_pair] = 1;
//     }

//     std::cout << "\n\n==============[NORMALIZATION: first parent]=============\n\n";
//     for(const auto& elem : first_parent_colour_occurrences)
//     {
//         std::cout << "\nColour: " << elem.first << ", count: " << elem.second;
//     }

//     std::cout << "\n\n==============[NORMALIZATION: second parent]=============\n\n";
//     for(const auto& elem : second_parent_colour_occurrences)
//     {
//         std::cout << "\nColour: " << elem.first << ", count: " << elem.second;
//     }

//     std::cout << "\n\n==============[NORMALIZATION: common colours]=============\n\n";
//     for(const auto& elem : common_colours)
//     {
//         std::cout << "\nColours pair: [" << get<0>(elem.first) << ", " << get<1>(elem.first) << "], count: " << elem.second;
//     }
// }

// void GeneticSolver::evolve_population(
//     std::vector<Solution>& parents, std::vector<Solution>& offsprings)
// {
//     std::size_t current_position {0uz};

//     while (current_position < parents.size()) {
//         std::swap(population.at(current_position), parents.at(current_position));

        current_position++;
    }

    const std::size_t remaining_places {model_ref.model_params.vertices - current_position};
    std::size_t offspring_position {0uz};

    while (offspring_position < offsprings.size()) {
        if (offspring_position >= remaining_places)
            return;

        std::swap(population.at(current_position), offsprings.at(offspring_position));
        current_position++;
        offspring_position++;
    }
}

void GeneticSolver::process_mutation()
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    std::uniform_int_distribution<std::size_t> vertex_int_distribution(0, model_ref.model_params.vertices - 1);
    std::uniform_int_distribution<std::size_t> colour_int_distribution(1, model_ref.model_params.max_degree);

    for (auto& solution : population) {
        const double probability {distribution(model_ref.rng)};

        if (probability < mutation_probability) {
            const std::size_t random_vertex {vertex_int_distribution(model_ref.rng)};
            const std::size_t random_colour {colour_int_distribution(model_ref.rng)};

            solution.graph.vertices.at(random_vertex).update_colour(random_colour);
        }
    }
}

Solution GeneticSolver::create_new_solution(Graph&& graph)
{
    Solution solution(std::move(graph));
    solution.fitness = fitness_evaluation(solution);

    return solution;
}

double GeneticSolver::fitness_evaluation(Solution& solution)
{
    fitness_evaluations++;
    return model_ref.evaluate_fitness(solution.graph);
}
