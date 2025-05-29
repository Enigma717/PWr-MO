#include "subpopulation.hpp"
#include "model.hpp"
#include "utility_operators.hpp"

#include <iostream>
#include <vector>
#include <map>
#include <ranges>

namespace
{
    std::size_t tournament_size {2uz};
    std::size_t candidates_pair_step {2uz};

    constexpr std::pair<std::size_t, std::size_t> get_second(std::pair<std::size_t, std::size_t> elem)
    {
        return std::make_pair(elem.second, elem.first);
    }
}

std::size_t  Subpopulation::fitness_evaluations {0uz};

Subpopulation::Subpopulation(std::size_t subpopulation_size, Model& model_ref)
: subpopulation_size {subpopulation_size},
  best_solution {nullptr},
  worst_solution {nullptr},
  model_ref {model_ref}
{
    individuals.resize(subpopulation_size);
    improving_offsprings.reserve(subpopulation_size);

    for (auto& solution : individuals)
        solution = create_new_solution(model_ref.solve_random());
}

std::size_t Subpopulation::get_ffe()
{
    return fitness_evaluations;
}

void Subpopulation::print_individuals() const
{
    for (std::size_t i {0uz}; i < individuals.size(); i++)
        std::cout << "Subpopulation size: " << subpopulation_size << " | Solution " << i
            << " (" << &individuals.at(i) << "): " << individuals.at(i) << "\n";
}

void Subpopulation::print_info() const
{
    std::cout << "Subpopulation size: [" << subpopulation_size << "]"
        << " | Is locked: [" << is_locked << "]"
        << " | Best solution: [" << best_solution << "], fitness: " << best_solution->fitness
        << " | Worst solution: [" << worst_solution << "], fitness: " << worst_solution->fitness
        << " | Average fitness: " << avg_fitness
        << " | Iterations done: " << iterations_done << "\n";
}

void Subpopulation::run_iteration()
{
    update_subpopulation_data();

    auto candidates {tournament_selection()};
    process_crossover(candidates);
    subsitute_subpopulation_with_offsprings();

    iterations_done++;
}

Solution Subpopulation::create_new_solution(Graph&& graph)
{
    Solution solution(std::move(graph));
    solution.fitness = fitness_evaluation(solution);

    return solution;
}

double Subpopulation::fitness_evaluation(Solution& solution)
{
    fitness_evaluations++;
    return model_ref.evaluate_fitness(solution.graph);
}

void Subpopulation::subsitute_subpopulation_with_offsprings()
{
    const auto cutoff_point {std::min(improving_offsprings.size(), subpopulation_size)};

    if (cutoff_point == 0uz)
        return;

    for (std::size_t i {0uz}; i < cutoff_point; i++)
        individuals.at(i) = std::move(improving_offsprings.at(i));
}

void Subpopulation::update_subpopulation_data()
{
    auto new_best_solution {&*std::min_element(individuals.begin(), individuals.end())};
    auto new_worst_solution {&*std::max_element(individuals.begin(), individuals.end())};
    double new_avg_fitness {
        static_cast<double>(std::accumulate(individuals.begin(), individuals.end(), 0.0)) /
        static_cast<double>(subpopulation_size)};

    if (iterations_done == 0) {
        best_solution = new_best_solution;
        worst_solution = new_worst_solution;
        avg_fitness = new_avg_fitness;

        return;
    }

    if (new_best_solution->fitness <= best_solution->fitness)
        best_solution = new_best_solution;

    if (new_worst_solution->fitness > worst_solution->fitness)
        worst_solution = new_worst_solution;

    avg_fitness = new_avg_fitness;
}

std::vector<Solution*> Subpopulation::tournament_selection()
{
    std::vector<Solution*> final_winners;
    final_winners.reserve(subpopulation_size);

    std::uniform_int_distribution<std::size_t> int_distribution(0, subpopulation_size - 1);

    for (std::size_t i {0uz}; i < subpopulation_size; i++) {
        std::vector<Solution*> tournament_group;
        tournament_group.reserve(tournament_size);

        for (std::size_t j {0uz}; j < tournament_size; j++) {
            const auto generated_index {int_distribution(model_ref.rng)};
            const auto candidate {&individuals.at(generated_index)};

            tournament_group.push_back(candidate);
        }

        const auto winner {*std::min_element(tournament_group.begin(), tournament_group.end())};
        final_winners.push_back(winner);
    }

    return final_winners;
}

void Subpopulation::process_crossover(std::vector<Solution*> candidates)
{
    improving_offsprings.clear();

    for (std::size_t i {0uz}; i < candidates.size(); i += candidates_pair_step) {
        Solution* first_parent {candidates.at(i)};
        Solution* second_parent {candidates.at(i + 1)};

        if (first_parent != second_parent)
            process_partition_crossover({first_parent, second_parent});
    }
}

void Subpopulation::process_partition_crossover(const std::vector<Solution*>& parents)
{
    auto& first_parent_graph {parents.at(0)->graph};
    auto& second_parent_graph {parents.at(1)->graph};
    const auto building_blocks {normalize_parent_colours(first_parent_graph, second_parent_graph)};

    for (const auto& block : building_blocks) {
        Graph new_offspring_graph {second_parent_graph};

        for (const auto index : block) {
            const auto first_parent_colour {first_parent_graph.vertices.at(index).get_colour()};
            new_offspring_graph.vertices.at(index).update_colour(first_parent_colour);
        }

        const auto offspring {create_new_solution(std::move(new_offspring_graph))};

        if (offspring.fitness <= parents[0]->fitness
            && offspring.fitness <= parents[1]->fitness
            && offspring.fitness <= best_solution->fitness
        )
            improving_offsprings.push_back(std::move(offspring));
    }
}

BuildingBlocks Subpopulation::normalize_parent_colours(
    Graph& first_parent_graph,
    Graph& second_parent_graph)
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
        auto matching_pairs {common_colours_rates | std::views::filter([&](auto& v) {
            return v.first.first == pair.second && free_colours.contains(v.first.second);
        })};

        if (matching_pairs.empty())
            continue;

        const auto best_match {*std::max_element(
            matching_pairs.begin(),
            matching_pairs.end(),
            [](const auto& p1, const auto& p2) {
            return p1.second < p2.second;
        })};

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

    std::vector<std::size_t> mismatches;

    for (std::size_t i {0uz}; i < second_parent_graph.vertices.size(); i++) {
        Vertex& first_vertex {first_parent_graph.vertices.at(i)};
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

        if (first_vertex.get_colour() != second_vertex.get_colour())
            mismatches.push_back(i);
    }

    std::vector<bool> mismatches_used(mismatches.size());
    std::vector<std::vector<std::size_t>> building_blocks;

    for (std::size_t i {0uz}; i < mismatches.size(); i++) {
        if (mismatches_used.at(i))
            continue;

        const auto mismatch {mismatches.at(i)};
        std::vector<std::size_t> building_block;
        building_block.push_back(mismatch);
        mismatches_used.at(i) = true;

        auto neighbours {first_parent_graph.vertices.at(mismatch).get_neighbours()};

        for (auto neighbour : neighbours) {
            const auto it {std::find_if(
                mismatches.begin(),
                mismatches.end(),
                [=](const auto index) {
                    return index == neighbour->get_id();
            })};

            auto index {it - mismatches.begin()};
            if (it != mismatches.end() && !mismatches_used.at(index)) {
                building_block.push_back(*it);
                mismatches_used.at(index) = true;
            }
        }

        building_blocks.push_back(building_block);
    }

    return building_blocks;
}