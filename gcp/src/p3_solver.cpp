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
    constexpr std::size_t subpopulation_base_size {1uz};
    constexpr std::size_t ltgomea_base_iterations {4uz};
    constexpr std::size_t iterations_limit {100'000uz};
    constexpr std::size_t ffe_limit {500'000uz};

    constexpr std::pair<std::size_t, std::size_t> get_second(std::pair<std::size_t, std::size_t> elem)
    {
        return std::make_pair(elem.second, elem.first);
    }
}

P3Solver::P3Solver(Model& model_ref) : model_ref{model_ref} {}

void P3Solver::print_info() const
{
    std::cout << "Pyramid levels: [" << pyramid_levels << "]"
        << " | Best solution: [" << best_solution << "], fitness: " << best_solution->fitness
        << " | Worst solution: [" << worst_solution << "], fitness: " << worst_solution->fitness
        << " | Average fitness: " << avg_fitness
        << " | Iterations done: " << total_iterations << "\n";
}

void P3Solver::solve()
{
    total_iterations = 0uz;

    create_new_level();

    // while (!is_optimum_reached
        // && total_iterations < iterations_limit
        // && Subpopulation::get_ffe() < ffe_limit
    // ) {
    for (size_t i {0uz}; i < 100; i++)
        next_iteration();
    // }

    std::cout << "\n\nFinal map (" << known_solutions.size() << "):";
    for (const auto& solution : known_solutions) {
        std::cout << "\nSolution [" << solution.first << "]: " << solution.second;
    }

    std::cout << "\n\nFinal pyramid (" << pyramid.size() << "):";
    for (std::size_t level {0uz}; level < pyramid_levels; level++) {
        const auto& pyramid_level {pyramid.at(level)};
        std::cout << "\nLevel number: " << level << " | Level size: " << pyramid_level.size();

        for (const auto& solution : pyramid_level)
            std::cout << "\nSolution [" << &solution << "]: " << solution.fitness;
    }

    // Solution* final_solution {best_solution};
//
    // std::cout << "\n\nFinal results:\n";
    // for (const auto& level : pyramid) {
        // for (const auto& individual : level)
            // if (level.best_solution->fitness < final_solution->fitness)
                // final_solution = level.best_solution;
//
        // level.print_info();
    // }

    // std::cout << "\nFinal solution: [" << final_solution << "] | " << *final_solution;
    std::cout << "\nUsed crossover type: [" << (crossover_type == 0u ? "optimal_mixing]\n" : "partition_crossover]\n");

    return;
}

void P3Solver::create_new_level()
{
    pyramid.push_back({});
    linkage_trees.push_back({model_ref.base_graph->vertices.size()});
    pyramid_levels++;
}

void P3Solver::add_solution_to_level(const Solution& solution, const std::size_t level)
{
    if (level >= pyramid_levels)
        create_new_level();

    auto& current_level {pyramid.at(level)};
    auto& current_linkage_tree {linkage_trees.at(level)};

    current_level.push_back(solution);
    current_linkage_tree.calculate_DSM(current_level);
    current_linkage_tree.create_clusters();

    // std::cout << "\nClusters: [\n";
    // for (const auto& cluster : linkage_trees.at(level).clusters)
        // std::cout << " [" << cluster << "]";
    // std::cout << " ]";
}

void P3Solver::next_iteration()
{
    Solution new_solution(create_new_solution(model_ref.solve_random()));

    if (pyramid.at(0).size() > 0)
        normalize_colours(pyramid.at(0).at(0).graph, new_solution.graph);

    if (!known_solutions.contains(std::hash<Solution>{}(new_solution))) {
        known_solutions.insert({std::hash<Solution>{}(new_solution), new_solution.fitness});
        add_solution_to_level(new_solution, 0);
    }
    // else {
        // std::cout << "\n[Y]Solution [" << new_solution << "] found in map";
    // }

    for (std::size_t current_level {0uz}; current_level < pyramid_levels; current_level++) {
        process_optimal_mixing(new_solution, current_level);

        if (!known_solutions.contains(std::hash<Solution>{}(new_solution))) {
            std::cout << "\n[XD] Solution [" << &new_solution << "] NOT found in map, fitness: " << new_solution.fitness;
            known_solutions.insert({std::hash<Solution>{}(new_solution), new_solution.fitness});
            add_solution_to_level(new_solution, current_level + 1);
        }
        // else {
            // std::cout << "\n[XD]Solution [" << new_solution << "] found in map";
        // }
    }

    total_iterations++;
}

Solution P3Solver::create_new_solution(Graph&& graph)
{
    Solution solution(std::move(graph));
    solution.fitness = fitness_evaluation(solution);

    return solution;
}

double P3Solver::fitness_evaluation(Solution& solution)
{
    fitness_evaluations++;
    return model_ref.evaluate_fitness(solution.graph);
}

void P3Solver::normalize_colours(Graph& first_parent_graph, Graph& second_parent_graph)
{
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
        const double rate {
            static_cast<double>(elem.second) /
            static_cast<double>(first_parent_colour_occurrences[elem.first.first])};
        common_colours_rates[elem.first] = rate;
    }

    std::vector<ColoursPair> colours_order(first_parent_colour_occurrences.size());
    std::transform(
        first_parent_colour_occurrences.begin(),
        first_parent_colour_occurrences.end(),
        colours_order.begin(),
        get_second);
    std::sort(colours_order.begin(), colours_order.end(), std::greater{});

    std::set<std::size_t> free_colours;
    for (const auto& colours : second_parent_colour_occurrences)
        free_colours.insert(colours.first);

    std::set<std::size_t> unmatched_colours;
    for (const auto& colours : first_parent_colour_occurrences)
        unmatched_colours.insert(colours.first);


    std::vector<ColoursPair> colours_to_use;
    colours_to_use.reserve(first_parent_colour_occurrences.size());

    for (const auto& pair : colours_order) {
        auto matching_pairs = common_colours_rates | std::views::filter([&](auto& v) {
            return v.first.first == pair.second && free_colours.contains(v.first.second);
        });

        if (matching_pairs.empty())
            continue;

        auto best_match = *std::max_element(
            matching_pairs.begin(),
            matching_pairs.end(),
            [](const auto& p1, const auto& p2) {
            return p1.second < p2.second;
        });

        free_colours.erase(best_match.first.second);
        unmatched_colours.erase(best_match.first.first);

        colours_to_use.push_back(std::make_pair(best_match.first.first, best_match.first.second));
    }

    if (!free_colours.empty() && !unmatched_colours.empty()) {
        for (const auto unmatched_colour : unmatched_colours) {
            auto first_free_colour {*free_colours.begin()};
            colours_to_use.push_back(std::make_pair(unmatched_colour, first_free_colour));
            free_colours.erase(first_free_colour);
        }
    }

    for (std::size_t i {0uz}; i < second_parent_graph.vertices.size(); i++) {
        Vertex& second_vertex {second_parent_graph.vertices.at(i)};
        const std::size_t second_vertex_colour {second_vertex.get_colour()};

        const auto new_colour_pair {std::find_if(
            colours_to_use.begin(),
            colours_to_use.end(),
            [=](const auto& pair) {
                return second_vertex_colour == pair.second;
        })};

        if (new_colour_pair != colours_to_use.end())
            second_vertex.update_colour(new_colour_pair->first);
    }
}

void P3Solver::process_optimal_mixing(Solution& solution, std::size_t current_level)
{
    const auto& clusters {linkage_trees.at(current_level).clusters};

    if (clusters.size() == 0)
        return;

    // std::cout << "\n========[NOWA ITERACJA]========\n";
    const auto& pyramid_level {pyramid.at(current_level)};

    std::uniform_int_distribution<std::size_t> int_distribution(0, pyramid_level.size() - 1);
    const auto generated_index {int_distribution(model_ref.rng)};

    // std::cout << "\nGenerated index: " << generated_index;
    Solution backup = solution;
    bool skip_cluster {false};

    // std::cout << "\nGOM clusters: [\n";
    for (const auto& cluster : clusters) {
        skip_cluster = false;
        // std::cout << "\n[" << cluster << "]\n";

        auto* donor {&pyramid_level.at(generated_index)};

        std::vector<bool> donors_tried(pyramid_level.size());
        donors_tried.at(generated_index) = true;

        // std::cout << "\n[BEFORE] Donors_tried: [";
        // for (const auto xd : donors_tried)
            // std::cout << " " << xd;
        // std::cout << "]\n";

        // print_individuals();
        // std::cout << "\n[INSIDE]Changed individual colors: [" << current_offspring.graph.colours << "]";
        // std::cout << "\n[INSIDE]Drawn donor colors: [" << donor->graph.colours << "]";
        while (model_ref.check_for_equality_in_cluster(solution, *donor, cluster)) {
            // std::cout << "\nChanged individual and donor are the same: ";

            if (std::find(donors_tried.begin(), donors_tried.end(), false) == donors_tried.end()) {
                // std::cout << "\nNo more unique donors for this individual";
                skip_cluster = true;
                break;
            }
            // std::cout << "\n[INSIDE]Changed individual: [" << &current_offspring << "]: " << current_offspring;
            // std::cout << "\n[INSIDE]Drawn donor: [" << donor << "]: " << *donor;

            const auto new_index {int_distribution(model_ref.rng)};
            donors_tried.at(new_index) = true;
            // std::cout << "\nNew index: " << new_index;
            // std::cout << "\nIndividuals size: " << individuals.size();
            donor = &pyramid_level.at(new_index);
            // std::cout << "\n[NEW]Changed individual: [" << &changed << "]: " << changed;
            // std::cout << "\n[NEW]Drawn donor: [" << donor << "]: " << *donor;

            // std::cout << "\n[AFTER] Donors_tried: [";
            // for (const auto xd : donors_tried)
            //     std::cout << " " << xd;
            // std::cout << "]\n";
        }

        if (skip_cluster)
            continue;

        // std::cout << "\nChanged individual: [" << &changed << "]: " << changed;
        // std::cout << "\nBackup: [" << &backup << "]: " << backup;
        // std::cout << "\nDrawn donor: [" << donor << "]: " << *donor;

        // std::cout << "\nColors from ind and donor for given cluster: ";
        // std::cout << "\nChanged: [ ";
        // for (const auto value : cluster) {
        //     std::cout << " " << *changed.graph.colours.at(value) << ",";
        // }
        // std::cout << " ]";

        // std::cout << "\nBackup: [ ";
        // for (const auto value : cluster) {
        //     std::cout << " " << *backup.graph.colours.at(value) << ",";
        // }
        // std::cout << " ]";

        // std::cout << "\nDonor: [ ";;
        // for (const auto value : cluster) {
        //     std::cout << " " << *donor->graph.colours.at(value) << ",";
        // }
        // std::cout << " ]";


        // std::cout << "\n[BEFORE] Donor: [" << donor << "]: " << *donor;
        // std::cout << "\n[BEFORE] Offpsring: [" << &solution << "]: " << solution;
        // std::cout << "\n[BEFORE] Backup: [" << &backup << "]: " << backup;

        // std::cout << "\nTrying to donate";
        for (const auto value : cluster) {
            solution.graph.vertices.at(value).colour = donor->graph.vertices.at(value).colour;
        }

        solution.fitness = fitness_evaluation(solution);

        // std::cout << "\n[MEANWHILE] Donor: [" << donor << "]: " << *donor;
        // std::cout << "\n[MEANWHILE] Offspring: [" << &solution << "]: " << solution;
        // std::cout << "\n[MEANWHILE] Backup: [" << &backup << "]: " << backup;

        if (solution.fitness > backup.fitness) {
            // std::cout << "\nWOW! Dubstep remix!\n";
            solution = backup;
        }
        else {
            backup = solution;
        }

        // std::cout << "\n[AFTER] Donor: [" << donor << "]: " << *donor;
        // std::cout << "\n[AFTER] Offspring: [" << &solution << "]: " << solution;
        // std::cout << "\n[AFTER] Backup: [" << &backup << "]: " << backup;
    }
    // std::cout << "\n]\n\n";


    // std::cout << "\nOffsprings at the end:";
    // for (const auto& offspring : offsprings) {
        // std::cout << "\nOffspring [ " << &offspring << "]: " << offspring;
    // }
    // std::cout << "\n";
}
