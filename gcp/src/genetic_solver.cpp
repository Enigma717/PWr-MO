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
#include <utility>

namespace
{
    constexpr std::size_t minimum_subgroups {1uz};
    constexpr std::size_t expected_winners {2uz};
    constexpr std::size_t dimension_divisor {4uz};
    constexpr std::size_t parents_pair_step {2uz};
    constexpr std::size_t print_info_threshold {100uz};
    constexpr std::size_t generation_limit {50'000uz};
    constexpr std::size_t ffe_limit {25'000uz};
    constexpr double random_initialization_probability {0.9};

    constexpr std::pair<std::size_t, std::size_t> get_second(std::pair<std::size_t, std::size_t> elem)
    {
        return std::make_pair(elem.second, elem.first);
    }
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
        << "\t||\tBest fitness: " << best_solution->fitness
        << " \t||\tWorst fitness: " << worst_solution->fitness
        << "\t||\tAverage fitness: " << avg_fitness
        << "\t||\tFFE: " << fitness_evaluations << "\n";

    return log.str();
}

bool GeneticSolver::check_reached_optimum()
{
    return best_solution->fitness == model_ref.model_params.optimum;
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
    auto new_best_solution {&*std::min_element(population.begin(), population.end())};
    auto new_worst_solution {&*std::max_element(population.begin(), population.end())};
    double new_avg_fitness {
        static_cast<double>(std::accumulate(population.begin(), population.end(), 0.0)) /
        static_cast<double>(population_size)};

    if (generation_number == 0) {
        best_solution = new_best_solution;
        worst_solution = new_worst_solution;
        avg_fitness = new_avg_fitness;
    }

    if (new_best_solution->fitness < best_solution->fitness) {
        std::swap(best_solution, new_best_solution);
        std::cout << print_generation_info();
    }

    if (new_worst_solution->fitness > worst_solution->fitness) {
        std::swap(worst_solution, new_worst_solution);
    }

    avg_fitness = new_avg_fitness;

    if (generation_number % print_info_threshold == 0 && generation_number != 0)
        std::cout << print_generation_info();
}

Solution& GeneticSolver::solve()
{
    initialize_population(PopulationType::random, population_size);
    // print_population();

    std::cout << "\n\n===========================\n\n";

    // Graph first_graph(10);
    // Graph second_graph(10);

    // first_graph.vertices.at(0).update_colour((11 - 1));
    // first_graph.vertices.at(1).update_colour((12 - 1));
    // first_graph.vertices.at(2).update_colour((13 - 1));
    // first_graph.vertices.at(3).update_colour((14 - 1));
    // first_graph.vertices.at(4).update_colour((14 - 1));
    // first_graph.vertices.at(5).update_colour((11 - 1));
    // first_graph.vertices.at(6).update_colour((12 - 1));
    // first_graph.vertices.at(7).update_colour((11 - 1));
    // first_graph.vertices.at(8).update_colour((14 - 1));
    // first_graph.vertices.at(9).update_colour((11 - 1));

    // first_graph.add_edge(0, 1);
    // first_graph.add_edge(0, 2);
    // first_graph.add_edge(0, 3);
    // first_graph.add_edge(1, 0);
    // first_graph.add_edge(1, 2);
    // first_graph.add_edge(1, 5);
    // first_graph.add_edge(2, 0);
    // first_graph.add_edge(2, 1);
    // first_graph.add_edge(2, 3);
    // first_graph.add_edge(2, 6);
    // first_graph.add_edge(3, 0);
    // first_graph.add_edge(3, 2);
    // first_graph.add_edge(3, 5);
    // first_graph.add_edge(4, 7);
    // first_graph.add_edge(4, 9);
    // first_graph.add_edge(5, 1);
    // first_graph.add_edge(5, 3);
    // first_graph.add_edge(5, 9);
    // first_graph.add_edge(6, 2);
    // first_graph.add_edge(6, 4);
    // first_graph.add_edge(7, 4);
    // first_graph.add_edge(7, 8);
    // first_graph.add_edge(8, 7);
    // first_graph.add_edge(8, 9);
    // first_graph.add_edge(9, 4);
    // first_graph.add_edge(9, 5);
    // first_graph.add_edge(9, 8);

    // second_graph.vertices.at(0).update_colour(2);
    // second_graph.vertices.at(1).update_colour(2);
    // second_graph.vertices.at(2).update_colour(1);
    // second_graph.vertices.at(3).update_colour(3);
    // second_graph.vertices.at(4).update_colour(3);
    // second_graph.vertices.at(5).update_colour(4);
    // second_graph.vertices.at(6).update_colour(3);
    // second_graph.vertices.at(7).update_colour(3);
    // second_graph.vertices.at(8).update_colour(1);
    // second_graph.vertices.at(9).update_colour(2);

    // second_graph.add_edge(0, 1);
    // second_graph.add_edge(0, 2);
    // second_graph.add_edge(0, 3);
    // second_graph.add_edge(1, 0);
    // second_graph.add_edge(1, 2);
    // second_graph.add_edge(1, 5);
    // second_graph.add_edge(2, 0);
    // second_graph.add_edge(2, 1);
    // second_graph.add_edge(2, 3);
    // second_graph.add_edge(2, 6);
    // second_graph.add_edge(3, 0);
    // second_graph.add_edge(3, 2);
    // second_graph.add_edge(3, 5);
    // second_graph.add_edge(4, 7);
    // second_graph.add_edge(4, 9);
    // second_graph.add_edge(5, 1);
    // second_graph.add_edge(5, 3);
    // second_graph.add_edge(5, 9);
    // second_graph.add_edge(6, 2);
    // second_graph.add_edge(6, 4);
    // second_graph.add_edge(7, 4);
    // second_graph.add_edge(7, 8);
    // second_graph.add_edge(8, 7);
    // second_graph.add_edge(8, 9);
    // second_graph.add_edge(9, 4);
    // second_graph.add_edge(9, 5);
    // second_graph.add_edge(9, 8);

    // Solution first(std::move(first_graph));
    // Solution second(std::move(second_graph));

    // std::swap(population.at(0), first);
    // std::swap(population.at(1), second);

    evaluate_population();

    while (!check_reached_ffe_limit() && !check_reached_gen_limit() && !check_reached_optimum()) {
        std::vector<Solution> parents {tournament_selection(subgroup_size)};
        std::vector<Solution> offsprings {crossover_parents(parents)};

        evolve_population(parents, offsprings);
        process_mutation();
        evaluate_population();

        generation_number++;
    }

    std::cout << print_generation_info();

    return *best_solution;
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
    const std::size_t subgroups_count {std::max(minimum_subgroups, population.size() / subgroup_size)};
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

std::vector<Solution> GeneticSolver::crossover_parents(std::vector<Solution>& parents)
{
    std::uniform_real_distribution<double> real_distribution(0.0, 1.0);

    std::vector<Solution> offsprings;
    offsprings.reserve(parents.size());

    for (std::size_t i {0uz}; i < parents.size(); i += parents_pair_step) {
        const double probability {real_distribution(model_ref.rng)};
        Solution first_parent {parents.at(i)};
        Solution second_parent {parents.at(i + 1)};
        Graph& first_parent_graph {first_parent.graph};
        Graph& second_parent_graph {second_parent.graph};

        if (probability < crossing_probability) {
            const auto new_offsprings {process_crossover(first_parent_graph, second_parent_graph)};

            if (offsprings.size()) {
                for (const auto& offspring : new_offsprings)
                    offsprings.push_back(offspring);
            }
            else {
                offsprings.push_back(first_parent);
                offsprings.push_back(second_parent);
            }
        }
        else {
            offsprings.push_back(first_parent);
            offsprings.push_back(second_parent);
        }
    }

    return offsprings;
}

std::vector<Solution> GeneticSolver::process_crossover(
    Graph& first_parent_graph, Graph& second_parent_graph)
{
    switch (crossover_type) {
        case CrossoverType::dpoint:
            return process_double_point_crossover(first_parent_graph, second_parent_graph);
        case CrossoverType::uniform:
            return process_uniform_crossover(first_parent_graph, second_parent_graph);
        case CrossoverType::partition:
            return process_partition_crossover(first_parent_graph, second_parent_graph);
    }

    return {};
}

std::vector<Solution> GeneticSolver::process_double_point_crossover(
        const Graph& first_parent_graph, const Graph& second_parent_graph)
{
    const std::size_t dimension {model_ref.model_params.vertices};
    const std::size_t dimension_quarter {dimension / dimension_divisor};

    std::uniform_int_distribution<std::size_t> int_distribution(
        dimension_quarter, dimension - dimension_quarter);

    std::size_t first_crossing_point {int_distribution(model_ref.rng)};
    std::size_t second_crossing_point {int_distribution(model_ref.rng)};

    if (first_crossing_point > second_crossing_point)
        std::swap(first_crossing_point, second_crossing_point);

    const Solution first_offspring {double_point_crossover(
        first_parent_graph, second_parent_graph, first_crossing_point, second_crossing_point)};
    const Solution second_offspring {double_point_crossover(
        second_parent_graph, first_parent_graph, first_crossing_point, second_crossing_point)};

    return {first_offspring, second_offspring};
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

std::vector<Solution> GeneticSolver::process_uniform_crossover(
        const Graph& first_parent_graph, const Graph& second_parent_graph)
{
    const Solution first_offspring {uniform_crossover(first_parent_graph, second_parent_graph)};
    const Solution second_offspring {uniform_crossover(second_parent_graph, first_parent_graph)};

    return {first_offspring, second_offspring};
}

Solution GeneticSolver::uniform_crossover(
    const Graph& first_parent_graph,
    const Graph& second_parent_graph)
{
    const std::size_t dimension {first_parent_graph.vertices.size()};
    Graph offspring_graph {first_parent_graph};

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (std::size_t i {0uz}; i < dimension; i++) {
        const double probability {distribution(model_ref.rng)};

        if (probability <= 0.5)
            offspring_graph.vertices.at(i).update_colour(first_parent_graph.vertices.at(i).get_colour());
        else
            offspring_graph.vertices.at(i).update_colour(second_parent_graph.vertices.at(i).get_colour());
    }

    return create_new_solution(std::move(offspring_graph));
}

std::vector<Solution> GeneticSolver::process_partition_crossover(
    Graph& first_parent_graph,
    Graph& second_parent_graph)
{
    const auto building_blocks {normalize_parent_colours(first_parent_graph, second_parent_graph)};

    std::vector<Solution> offsprings;
    offsprings.reserve(building_blocks.size());

    // std::cout << "\n-------------[FINAL BUILDING BLOCKS]---------------\n";
    for (const auto& block : building_blocks) {
        // std::cout << "\n\nBlock in partition: [" << block << "]\n";
        Graph new_offspring_graph {second_parent_graph};

        // std::cout << "\nNew offspring graph before change: " << new_offspring_graph;
        for (const auto index : block) {
            const auto first_parent_colour {first_parent_graph.vertices.at(index).get_colour()};
            // std::cout << "\nUpdating offspring at vertex: " << index + 1 << " into colour: " << std::hex << first_parent_colour;
            new_offspring_graph.vertices.at(index).update_colour(first_parent_colour);
        }
        // std::cout << "\nNew offspring graph after change: " << new_offspring_graph;
        offsprings.push_back(create_new_solution(std::move(new_offspring_graph)));
        // std::cout << "\nNew offspring graph after pushing: " << offsprings.back();
    }

    // std::cout << "\n\n\nFINAL OFFSPRINGS: [";
    // for (const auto& solution : offsprings) {
        // std::cout << "\nOffspring: [" << solution << "]";
        // std::cout << "\nFIRST VERTEX COLOUR: " << solution.graph.vertices.at(0).get_colour();
    // }
    // std::cout << "\n]\n\n";

    return offsprings;
}

 BuildingBlocks GeneticSolver::normalize_parent_colours(
    Graph& first_parent_graph,
    Graph& second_parent_graph)
{
    // std::cout << "\n\n==============[NORMALIZATION]=============\n\n";
    // std::cout << "\nFirst parent (" << &first_parent_graph << "): " << first_parent_graph;
    // std::cout << "\nSecond parent (" << &second_parent_graph << "): " << second_parent_graph;

    using ColoursMap = std::map<std::size_t, std::size_t>;
    using ColoursPair = std::pair<std::size_t, std::size_t>;

    ColoursMap first_parent_colour_occurrences;
    ColoursMap second_parent_colour_occurrences;
    std::map<ColoursPair, std::size_t> common_colours_counts;
    std::map<ColoursPair, double> common_colours_rates;

    for (std::size_t i {0uz}; i < first_parent_graph.vertices.size(); i++) {
        std::size_t first_parent_colour {*first_parent_graph.colours.at(i)};
        std::size_t second_parent_colour {*second_parent_graph.colours.at(i)};
        ColoursPair colours_pair {std::make_tuple(first_parent_colour, second_parent_colour)};

        auto first_it(first_parent_colour_occurrences.find(first_parent_colour));
        if (first_it != first_parent_colour_occurrences.end())
           first_it->second++;
        else
           first_parent_colour_occurrences[first_parent_colour] = 1;

        auto second_it(second_parent_colour_occurrences.find(second_parent_colour));
        if (second_it != second_parent_colour_occurrences.end())
           second_it->second++;
         else
           second_parent_colour_occurrences[second_parent_colour] = 1;

        auto pair_it(common_colours_counts.find(colours_pair));
        if (pair_it != common_colours_counts.end())
           pair_it->second++;
        else
           common_colours_counts[colours_pair] = 1;
    }

    for (const auto& elem : common_colours_counts) {
        const double rate {static_cast<double>(elem.second) / static_cast<double>(first_parent_colour_occurrences[elem.first.first])};
        // std::cout << "\nRate: " << rate;
        common_colours_rates[elem.first] = rate;
    }

    std::vector<ColoursPair> colours_order(first_parent_colour_occurrences.size());
    std::transform(
        first_parent_colour_occurrences.begin(),
        first_parent_colour_occurrences.end(),
        colours_order.begin(),
        get_second);
    std::sort(colours_order.begin(), colours_order.end(), std::greater{});


    // std::cout << "\n\n--------------[NORMALIZATION: first parent]--------------\n\n";
    // for(const auto& elem : first_parent_colour_occurrences)
        // std::cout << "\nColour: " << std::hex << elem.first << ", count: " << std::dec << elem.second;
//
    // std::cout << "\n\n--------------[NORMALIZATION: second parent]--------------\n\n";
    // for(const auto& elem : second_parent_colour_occurrences)
        // std::cout << "\nColour: " << std::hex << elem.first << ", count: " << std::dec << elem.second;
//
    // std::cout << "\n\n--------------[NORMALIZATION: common counts]--------------\n\n";
    // for(const auto& elem : common_colours_counts)
        // std::cout << "\nColours counts pairs: [" << std::hex << elem.first.first << ", " << std::dec << elem.first.second << "], count: " << elem.second;
//
    // std::cout << "\n\n--------------[NORMALIZATION: common rates]--------------\n\n";
    // for(const auto& elem : common_colours_rates)
        // std::cout << "\nColours rates pairs: [" << std::hex << elem.first.first << ", " << std::dec << elem.first.second << "], rate: " << elem.second;
//
    // std::cout << "\n\n--------------[NORMALIZATION: colours order]--------------\n\n";
    // for (const auto& pair : colours_order)
        // std::cout << "\nColours orders pairs: [" << pair.first << ", " << std::hex << pair.second << "]";

    std::set<std::size_t> free_colours;
    for (const auto& colours : second_parent_colour_occurrences) {
        free_colours.insert(colours.first);
    }

    std::set<std::size_t> unmatched_colours;
    for (const auto& colours : first_parent_colour_occurrences) {
        unmatched_colours.insert(colours.first);
    }


    std::vector<ColoursPair> colours_to_use;
    colours_to_use.reserve(first_parent_colour_occurrences.size());

    std::size_t iteration {0uz};

    for (const auto& pair : colours_order) {
        // std::cout << "\nIteration: " << iteration;
        // std::cout << "\nFree colours in for loop start: [" << free_colours << "]";
        // std::cout << "\nUnmatched colours in for loop start: [" << unmatched_colours << "]";

        //FILTER THE ELEMENTS
        auto matching_pairs = common_colours_rates | std::views::filter([&](auto& v) {
            // std::cout << "\n\nINSIDE FILTER V.FIRST.FIRST: " << v.first.first;
            // std::cout << "\nINSIDE FILTER PAIR.SECOND: " << pair.second;
            // std::cout << "\nINSIDE FILTER V.FIRST SECOND: " << v.first.second;
            // std::cout << "\nINSIDE FILTER CONTAINS: " << free_colours.contains(v.first.second) << "\n";
            return v.first.first == pair.second && free_colours.contains(v.first.second);
        });

        if (matching_pairs.empty()) {
            // std::cout << "\nMATHICNG PAIRS EMPTY1 :DDD\n";
//
            // std::cout << "\nFree colours: [" << free_colours << "]";
            // std::cout << "\nUnmatched colours: [" << unmatched_colours << "]";
//
            // std::cout << "\nMATHICNG PAIRS EMPTY2 :DDD\n";
            continue;
        }

        // for (auto m : matching_pairs)
            // std::cout << "\nFound match: [" << m.first.first << ", " << m.first.second << "], rate: " << m.second;

        auto best_match = *std::max_element(
            matching_pairs.begin(),
            matching_pairs.end(),
            [](const auto& p1, const auto& p2) {
            return p1.second < p2.second;
        });

        free_colours.erase(best_match.first.second);
        unmatched_colours.erase(best_match.first.first);

        // std::cout << "\nBest match: [" << best_match.first.first << ", " << best_match.first.second << "], rate: " << best_match.second;

        colours_to_use.push_back(std::make_pair(best_match.first.first, best_match.first.second));
        iteration++;

        // std::cout << "\nFree colours in for loop end: [" << free_colours << "]";
        // std::cout << "\nUnmatched colours in for loop end: [" << unmatched_colours << "]";
    }

    // std::cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SIEMA ENUI:D \n";

    // std::cout << "\nFree colours: [" << free_colours << "]";
    // std::cout << "\nUnmatched colours: [" << unmatched_colours << "]";

    if (!free_colours.empty() && !unmatched_colours.empty()) {
        // std::cout << "\nŁooo tego1";
        // std::cout << "\nNot empty";

        for (const auto unmatched_colour : unmatched_colours) {
            // std::cout << "\nŁooo tego2";
            auto first_free_colour {*free_colours.begin()};
            // std::cout << "\nŁooo tego3";
            colours_to_use.push_back(std::make_pair(unmatched_colour, first_free_colour));
            // std::cout << "\nŁooo tego4: " << first_free_colour;
            free_colours.erase(first_free_colour);
            // std::cout << "\nŁooo tego5: " << unmatched_colour;
            // unmatched_colours.erase(unmatched_colour);
            // std::cout << "\nŁooo tego6";
        }
        // std::cout << "\nŁooo tego7";
    }

    // std::cout << "\nŁooo tego8";
    // std::cout << "\nFree colours: [" << free_colours << "]";
    // std::cout << "\nUnmatched colours: [" << unmatched_colours << "]";


    // for (const auto& best_pair : colours_to_use)
        // std::cout << "\nMatching pair: [" << best_pair.first << ", " << best_pair.second << "]";

    // std::cout << "\nFirst parent before normalization: " << first_parent_graph;
    // std::cout << "\nSecond parent before normalization: " << second_parent_graph;

    std::vector<std::size_t> mismatches;

    for (std::size_t i {0uz}; i < second_parent_graph.vertices.size(); i++) {
        Vertex& first_vertex {first_parent_graph.vertices.at(i)};
        Vertex& second_vertex {second_parent_graph.vertices.at(i)};
        const std::size_t second_vertex_colour {second_vertex.get_colour()};

        const auto new_colour_pair {std::find_if(colours_to_use.begin(), colours_to_use.end(), [=](const auto& pair) {return second_vertex_colour == pair.second;})};

        if (new_colour_pair != colours_to_use.end())
            second_vertex.update_colour(new_colour_pair->first);

        if (first_vertex.get_colour() != second_vertex.get_colour()) {
            // std::cout << "\nColours mismatch at vertex [" << i + 1 << "]: first parent - " << first_vertex.get_colour() << " | second parent - " << second_vertex.get_colour();
            mismatches.push_back(i);
        }
    }

    // std::cout << "\nFirst parent after normalization: " << first_parent_graph;
    // std::cout << "\nSecond parent after normalization: " << second_parent_graph;

    std::vector<bool> mismatches_used(mismatches.size());
    std::vector<std::vector<std::size_t>> building_blocks;

    for (std::size_t i {0uz}; i < mismatches.size(); i++) {
        if (mismatches_used.at(i))
            continue;

        const auto mismatch {mismatches.at(i)};
        // std::cout << "\nCurrent mismatch: " << mismatch;
        // std::cout << "\nMismatches: [" << mismatches << "]";
        // std::cout << "\nMismatches_used: [";
        // for (const auto used : mismatches_used)
            // std::cout << used << ", ";
        // std::cout << "]";
        std::vector<std::size_t> building_block;
        building_block.push_back(mismatch);
        mismatches_used.at(i) = true;

        auto neighbours {first_parent_graph.vertices.at(mismatch).get_neighbours()};

        // std::cout << "\n AFTER Mismatches: [" << mismatches << "]";
        // std::cout << "\n AFTER Mismatches_used: [";
        // for (const auto used : mismatches_used)
            // std::cout << used << ", ";
        // std::cout << "]";

        // std::cout << "\nNeighbours: " << neighbours;

        for (auto neighbour : neighbours) {
            // std::cout << "\nNeighbour get id: " << neighbour->get_id();
            const auto it {std::find_if(mismatches.begin(), mismatches.end(), [=](const auto index) {return index == neighbour->get_id();})};
            if (it != mismatches.end()) {
                // std::cout << "\nMAMY TO: " << *neighbour;
                std::size_t index {it - mismatches.begin()};
                building_block.push_back(*it);
                mismatches_used.at(index) = true;
            }
        }

        // std::cout << "\nBuilding block [" << building_block << "]\n";
        building_blocks.push_back(building_block);
    }

    // std::cout << "\n\nFINAL BUILDING BLOCKS: [";
    // for (const auto& block : building_blocks) {
        // std::cout << "\nBuilding block: [" << block << "]";
    // }

    // std::cout << "\n]\n\n";

    return building_blocks;
}

void GeneticSolver::evolve_population(
    std::vector<Solution>& parents, std::vector<Solution>& offsprings)
{
    std::size_t current_position {0uz};

    while (current_position < parents.size()) {
        std::swap(population.at(current_position), parents.at(current_position));
        current_position++;
    }

    const std::size_t remaining_places {model_ref.model_params.vertices - current_position};
    std::size_t offspring_position {0uz};

    while (offspring_position < offsprings.size() && current_position < population_size) {
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

    for (auto& solution : population) {
        const double probability {distribution(model_ref.rng)};

        if (probability < mutation_probability) {
            model_ref.mutate_random_vertex(solution.graph);
            solution.fitness = model_ref.evaluate_fitness(solution.graph);
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
