#include "../include/solver.h"
#include "../include/model.h"
#include "../include/graph.h"
#include "../include/vertex.h"
#include "../include/utility_operators.h"

#include <algorithm>
#include <iostream>
#include <vector>

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
