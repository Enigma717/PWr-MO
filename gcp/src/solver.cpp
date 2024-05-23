#include "../include/solver.h"
#include "../include/model.h"
#include "../include/graph.h"
#include "../include/vertex.h"
#include "../include/utility_operators.h"

#include <algorithm>
#include <iostream>

Solver::Solver(Model& model_ref) : model_ref{model_ref} {}

void Solver::random_solution() const
{
    for (auto& vertex : model_ref.graph->vertices)
        draw_colour(vertex);
}

void Solver::draw_colour(Vertex& vertex) const
{
    bool correctly_coloured {false};
    const auto& neighbours {vertex.get_neighbours()};

    do {
        std::uniform_int_distribution<std::size_t> int_distribution(1, model_ref.model_params.vertices);
        std::size_t drawn_colour {int_distribution(model_ref.rng)};

        auto match_colour {
            [=](const Vertex* const vertex) {
                return vertex->get_colour() == drawn_colour;
            }};

        auto it {std::find_if(
            neighbours.begin(),
            neighbours.end(),
            match_colour
        )};

        if (it == neighbours.end()) {
            vertex.update_colour(drawn_colour);
            correctly_coloured = true;
        }
    } while (!correctly_coloured);
}
