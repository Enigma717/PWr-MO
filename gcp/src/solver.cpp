#include "solver.hpp"
#include "model.hpp"
#include "graph.hpp"
#include "vertex.hpp"
#include "utility_operators.hpp"

#include <algorithm>
#include <iostream>
#include <chrono>
#include <vector>
#include <sstream>

namespace
{
    constexpr std::size_t maximum_sa_iterations {150'000uz};
}

Solver::Solver(Model& model_ref) : model_ref{model_ref} {}

void Solver::random_solution(Graph& graph) const
{
    for (auto& vertex : graph.vertices)
        colour_vertex_randomly(vertex);
}

void Solver::greedy_solution(Graph& graph) const
{
    for (auto& vertex : graph.vertices)
        colour_vertex_greedily(vertex);
}

void Solver::simulated_annealing_solution(Graph& graph) const
{
    std::ofstream results_file;
    std::stringstream results_path;
    results_path << "./csv/results/sa/sa_results_" << model_ref.model_params.instance_name << ".csv";
    results_file.open(results_path.str(), std::ios_base::app);

    std::ofstream plot_file;
    std::stringstream plot_data_path;
    plot_data_path << "./csv/results/sa/sa_plot_" << model_ref.model_params.instance_name << ".csv";
    plot_file.open(plot_data_path.str());

    if(!results_file.is_open()) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    if(!plot_file.is_open()) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    plot_file << "it; best; temp\n";

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    double current_temperature {sa_initial_temperature};

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    random_solution(graph);
    Graph best_graph {graph};
    Solution best_solution(std::move(best_graph));
    best_solution.fitness = model_ref.evaluate_fitness(best_solution.graph);

    std::size_t iteration {0};
    while (!check_reached_iteration_limit(iteration) && !check_reached_optimum(best_solution.fitness)) {
        Graph neighbour_graph {best_solution.graph};
        model_ref.mutate_random_vertex(neighbour_graph);
        Solution neighbour(std::move(neighbour_graph));
        neighbour.fitness = model_ref.evaluate_fitness(neighbour.graph);


        const auto fitness_diff {static_cast<int>(best_solution.fitness - neighbour.fitness)};
        const double probability {distribution(model_ref.rng)};
        const double exp {std::exp(fitness_diff / current_temperature)};

        if (fitness_diff >= 0) {
            std::cout << "|-> Iteration: " << iteration
                << "\t||\tBest fitness: " << best_solution.fitness
                << "\t||\tTemperature: " << current_temperature << "\n";

            std::swap(best_solution, neighbour);
        }
        else {
            if (probability < exp)
                std::swap(best_solution, neighbour);
        }

        if (iteration % 10 == 0) {
            plot_file << iteration << "; "
                << best_solution.fitness << "; "
                << current_temperature << "\n";
        }

        current_temperature *= sa_cooling_rate;
        iteration++;
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    results_file << sa_initial_temperature << "; "
        << sa_cooling_rate << "; "
        << best_solution.fitness << "; "
        << elapsed_time << "\n";

    plot_file << iteration << "; "
        << best_solution.fitness << "; "
        << current_temperature << "\n";

    std::swap(graph, best_solution.graph);

    plot_file.close();
    results_file.close();
}

bool Solver::check_reached_optimum(const std::size_t fitness) const
{
    return fitness == model_ref.model_params.optimum;
}

bool Solver::check_reached_iteration_limit(const std::size_t iteration) const
{
    return iteration >= maximum_sa_iterations;
}

void Solver::colour_vertex_randomly(Vertex& vertex) const
{
    bool correctly_coloured {false};
    std::size_t failed_retries {0uz};
    const auto forbidden_colours {model_ref.get_forbidden_colours(vertex)};

    do {
        std::uniform_int_distribution<std::size_t> int_distribution(1, model_ref.model_params.max_degree);
        std::size_t drawn_colour {int_distribution(model_ref.rng)};

        auto it {std::find_if(
            forbidden_colours.begin(),
            forbidden_colours.end(),
            [=](const std::size_t frobidden_colour) {
                return frobidden_colour == drawn_colour;
            }
        )};

        if (it == forbidden_colours.end()) {
            vertex.update_colour(drawn_colour);
            correctly_coloured = true;
        }
        else {
            failed_retries++;
        }

        if (failed_retries > (2 * model_ref.model_params.max_degree)) {
            vertex.update_colour(model_ref.model_params.max_degree + 1);
            correctly_coloured = true;
        }
    } while (!correctly_coloured);
}

void Solver::colour_vertex_greedily(Vertex& vertex) const
{
    const auto forbidden_colours {model_ref.get_forbidden_colours(vertex)};
    const std::size_t available_colour {model_ref.find_available_colour(forbidden_colours)};

    vertex.update_colour(available_colour);
}
