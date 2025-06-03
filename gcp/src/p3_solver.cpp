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
    constexpr std::size_t subpopulation_base_size {1uz};
    constexpr std::size_t ltgomea_base_iterations {4uz};
    constexpr std::size_t iterations_limit {100'000uz};
    constexpr std::size_t ffe_limit {1'000'000uz};

    constexpr std::pair<std::size_t, std::size_t> get_second(std::pair<std::size_t, std::size_t> elem)
    {
        return std::make_pair(elem.second, elem.first);
    }
}

P3Solver::P3Solver(Model& model_ref) : model_ref{model_ref} {}

void P3Solver::print_info() const
{
    std::cout << "Pyramid levels: " << pyramid_levels
        << " | Best solution fitness: " << best_solution.fitness
        << " | Average fitness: " << avg_fitness
        << " | Deviation: " << avg_deviation
        << " | Iterations: " << total_iterations
        << " | FFE: " << fitness_evaluations << "\n";
}

void P3Solver::solve()
{
    std::string dir_type {crossover_type == CrossoverType::optimal_mixing ? "p3" : "p3px"};

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

    results_file << "levels; best; avg; dev; iterations; ffe; time\n";
    plot_file << "iteration; levels; best; avg; dev\n";

    total_iterations = 0uz;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    create_new_level();

    while (best_solution.fitness != model_ref.model_params.optimum
        && total_iterations < iterations_limit
        && fitness_evaluations < ffe_limit
    ) {
        if (crossover_type == CrossoverType::optimal_mixing)
            next_om_iteration(plot_file);
        else if (crossover_type == CrossoverType::partition)
            next_px_iteration(plot_file);

        calculate_solver_info();
        print_info();
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    results_file << pyramid_levels << "; "
        << best_solution.fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "; "
        << total_iterations << "; "
        << fitness_evaluations << "; "
        << elapsed_time << "\n";

    plot_file << total_iterations << "; "
        << pyramid_levels << "; "
        << best_solution.fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "\n";

    plot_file.close();
    results_file.close();

    std::cout << "\n\nFinal solution: [" << &best_solution << "] | " << best_solution;
    std::cout << "\nIterations: " << total_iterations << " | fitness evaluations: " << fitness_evaluations;
    std::cout << "\nUsed crossover type: [" << (crossover_type == CrossoverType::optimal_mixing ? "optimal_mixing]\n" : "partition_crossover]\n");

    return;
}

void P3Solver::create_new_level()
{
    pyramid.push_back({});
    linkage_trees.push_back({model_ref.base_graph->vertices.size()});
    pyramid_levels++;
}

void P3Solver::calculate_solver_info()
{
    avg_fitness = 0.0;

    for (const auto& level : pyramid) {
        double level_avg_fitness {
            static_cast<double>(std::accumulate(level.begin(), level.end(), 0.0)) /
            static_cast<double>(level.size())};

        avg_fitness += level_avg_fitness;

        double sum {0.0};

        for (const auto& solution : level)
            sum += std::pow((solution.fitness - level_avg_fitness), 2);

        avg_deviation += sqrt(sum / level.size());
    }

    avg_fitness /= pyramid_levels;
    avg_deviation /= pyramid_levels;
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

void P3Solver::add_solution_to_level(Solution& solution, const std::size_t level)
{
    if (level >= pyramid_levels)
        create_new_level();

    if (solution.fitness < best_solution.fitness)
        best_solution = solution;

    auto& current_level {pyramid.at(level)};
    auto& current_linkage_tree {linkage_trees.at(level)};

    if (current_level.size() > 0)
        normalize_colours(current_level.at(0).graph, solution.graph);

    current_level.push_back(solution);
    current_linkage_tree.calculate_DSM(current_level);
    current_linkage_tree.create_clusters();
}

void P3Solver::normalize_colours(Graph& first_graph, Graph& second_graph)
{
    using ColoursMap = std::map<std::size_t, std::size_t>;
    using ColoursPair = std::pair<std::size_t, std::size_t>;

    ColoursMap first_parent_colour_occurrences;
    ColoursMap second_parent_colour_occurrences;
    std::map<ColoursPair, std::size_t> common_colours_counts;
    std::map<ColoursPair, double> common_colours_rates;

    for (std::size_t i {0uz}; i < first_graph.vertices.size(); i++) {
        std::size_t first_parent_colour {*first_graph.colours.at(i)};
        std::size_t second_parent_colour {*second_graph.colours.at(i)};
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

    for (std::size_t i {0uz}; i < second_graph.vertices.size(); i++) {
        Vertex& second_vertex {second_graph.vertices.at(i)};
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

void P3Solver::apply_hill_climber(Solution& solution)
{
    bool fitness_improved {false};
    const std::size_t vertices_count {solution.graph.vertices.size()};

    std::uniform_int_distribution<std::size_t> colour_distribution(1, model_ref.model_params.max_degree);

    std::unordered_map<std::size_t, std::size_t> processed_solutions;
    std::vector<std::size_t> walk_order(vertices_count);

    std::iota(walk_order.begin(), walk_order.end(), 0);

    do {
        fitness_improved = false;
        std::shuffle(walk_order.begin(), walk_order.end(), model_ref.rng);

        Solution backup {solution};

        for (const auto index : walk_order) {
            solution.graph.vertices.at(index).colour = colour_distribution(model_ref.rng);
            solution.fitness = fitness_evaluation(solution);

            const auto solution_hash {std::hash<Solution>{}(solution)};
            if (!processed_solutions.contains(solution_hash)) {
                processed_solutions.insert({solution_hash, solution.fitness});

                if (solution.fitness < backup.fitness)
                    fitness_improved = true;
                else
                    solution.graph.vertices.at(index).colour = backup.graph.vertices.at(index).colour;
            }
            else
                solution.graph.vertices.at(index).colour = backup.graph.vertices.at(index).colour;
        }
    } while (fitness_improved);
}

void P3Solver::next_om_iteration(std::ofstream& plot_file)
{
    Solution new_solution(create_new_solution(model_ref.solve_random()));
    apply_hill_climber(new_solution);

    if (total_iterations == 0)
        best_solution = new_solution;

    auto solution_hash {std::hash<Solution>{}(new_solution)};

    if (!known_solutions.contains(solution_hash)) {
        known_solutions.insert({solution_hash, new_solution.fitness});
        add_solution_to_level(new_solution, 0);
    }

    for (std::size_t current_level {0uz}; current_level < pyramid_levels; current_level++) {
        const auto previous_fitness {new_solution.fitness};
        process_optimal_mixing(new_solution, current_level);

        if (new_solution.fitness > previous_fitness)
            continue;

        const auto best_solution_in_level {
            &*std::min_element(pyramid.at(current_level).begin(), pyramid.at(current_level).end())};

        if (new_solution.fitness > best_solution_in_level->fitness)
            continue;

        solution_hash = std::hash<Solution>{}(new_solution);

        if (!known_solutions.contains(solution_hash)) {
            known_solutions.insert({solution_hash, new_solution.fitness});
            add_solution_to_level(new_solution, current_level + 1);
        }
    }

    plot_file << total_iterations << "; "
        << pyramid_levels << "; "
        << best_solution.fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "\n";

    total_iterations++;
}

void P3Solver::process_optimal_mixing(Solution& solution, std::size_t current_level)
{
    const auto& clusters {linkage_trees.at(current_level).clusters};
    const auto& pyramid_level {pyramid.at(current_level)};

    if (clusters.size() == 0)
        return;

    std::uniform_int_distribution<std::size_t> int_distribution(0, pyramid_level.size() - 1);
    const auto generated_index {int_distribution(model_ref.rng)};
    Solution backup {solution};
    bool skip_cluster {false};

    for (const auto& cluster : clusters) {
        skip_cluster = false;

        auto* donor {&pyramid_level.at(generated_index)};

        std::vector<bool> donors_tried(pyramid_level.size());
        donors_tried.at(generated_index) = true;

        while (model_ref.check_for_equality_in_cluster(solution, *donor, cluster)) {
            if (std::find(donors_tried.begin(), donors_tried.end(), false) == donors_tried.end()) {
                skip_cluster = true;
                break;
            }

            const auto new_index {int_distribution(model_ref.rng)};
            donors_tried.at(new_index) = true;
            donor = &pyramid_level.at(new_index);
        }

        if (skip_cluster)
            continue;

        for (const auto value : cluster)
            solution.graph.vertices.at(value).colour = donor->graph.vertices.at(value).colour;

        solution.fitness = fitness_evaluation(solution);

        if (solution.fitness > backup.fitness)
            solution = backup;
        else
            backup = solution;
    }
}

void P3Solver::next_px_iteration(std::ofstream& plot_file)
{
    Solution new_solution(create_new_solution(model_ref.solve_random()));
    apply_hill_climber(new_solution);

    if (total_iterations == 0)
        best_solution = new_solution;

    for (std::size_t current_level {0uz}; current_level < pyramid_levels; current_level++) {
        auto& pyramid_level {pyramid.at(current_level)};

        if (pyramid_level.size() <= 1) {
            const auto solution_hash {std::hash<Solution>{}(new_solution)};

            if (!known_solutions.contains(solution_hash)) {
                known_solutions.insert({solution_hash, new_solution.fitness});
                add_solution_to_level(new_solution, current_level);
            }

            plot_file << total_iterations << "; "
                << pyramid_levels << "; "
                << best_solution.fitness << "; "
                << avg_fitness << "; "
                << avg_deviation << "\n";

            total_iterations++;

            return;
        }

        std::uniform_int_distribution<std::size_t> partner_distribution(0, pyramid_level.size() - 2);
        auto* partner {&pyramid_level.at(partner_distribution(model_ref.rng))};

        const auto offspring {process_partition_crossover(new_solution, *partner)};

        if (offspring.fitness <= new_solution.fitness)
            new_solution = offspring;

        const auto best_solution_in_level {
            &*std::min_element(pyramid.at(current_level).begin(), pyramid.at(current_level).end())};

        if (new_solution.fitness > best_solution_in_level->fitness)
            continue;

        const auto solution_hash {std::hash<Solution>{}(new_solution)};

        if (!known_solutions.contains(solution_hash)) {
            known_solutions.insert({solution_hash, new_solution.fitness});
            add_solution_to_level(new_solution, current_level + 1);
        }
    }

    plot_file << total_iterations << "; "
        << pyramid_levels << "; "
        << best_solution.fitness << "; "
        << avg_fitness << "; "
        << avg_deviation << "\n";

    total_iterations++;
}

Solution P3Solver::process_partition_crossover(Solution& solution, Solution& partner)
{
    auto& solution_graph {solution.graph};
    auto& partner_graph {partner.graph};
    normalize_colours(partner_graph, solution_graph);

    const auto building_blocks {obtain_building_blocks(partner_graph, solution_graph)};

    std::vector<Solution> offsprings;
    offsprings.reserve(building_blocks.size());

    for (const auto& block : building_blocks) {
        Graph new_offspring_graph {partner_graph};

        for (const auto index : block) {
            const auto first_parent_colour {solution_graph.vertices.at(index).get_colour()};
            new_offspring_graph.vertices.at(index).update_colour(first_parent_colour);
        }

        offsprings.push_back(std::move(create_new_solution(std::move(new_offspring_graph))));
    }

    return offsprings.size() > 0 ? *std::min_element(offsprings.begin(), offsprings.end()) : solution;
}

BuildingBlocks P3Solver::obtain_building_blocks(
    const Graph& solution_graph,
    const Graph& partner_graph)
{
    std::vector<std::size_t> mismatches;

    for (std::size_t i {0uz}; i < partner_graph.vertices.size(); i++) {
        const Vertex& first_vertex {solution_graph.vertices.at(i)};
        const Vertex& second_vertex {partner_graph.vertices.at(i)};

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

        auto neighbours {solution_graph.vertices.at(mismatch).get_neighbours()};

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
