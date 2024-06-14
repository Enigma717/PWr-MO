#include "../include/solver.h"
#include "../include/model.h"
#include "../include/graph.h"
#include "../include/vertex.h"
#include "../include/utility_operators.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace
{
    constexpr double initial_sa_temperature {100'000.0};
    constexpr double cooling_rate {0.999};
    constexpr std::size_t maximum_sa_iterations {250'000uz};
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
    double current_temperature {initial_sa_temperature};

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


        const int fitness_diff {best_solution.fitness - neighbour.fitness};
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

        current_temperature *= cooling_rate;
        iteration++;
    }

    std::swap(graph, best_solution.graph);
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
    } while (!correctly_coloured);
}

void Solver::colour_vertex_greedily(Vertex& vertex) const
{
    const auto forbidden_colours {model_ref.get_forbidden_colours(vertex)};
    const std::size_t available_colour {model_ref.find_available_colour(forbidden_colours)};

    vertex.update_colour(available_colour);
}
