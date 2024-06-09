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
        colour_randomly(vertex);
}

void Solver::greedy_solution(Graph& graph) const
{
    for (auto& vertex : graph.vertices)
        vertex.update_colour(find_first_avail_colour(vertex));
}

void Solver::colour_randomly(Vertex& vertex) const
{
    bool correctly_coloured {false};
    const auto& neighbours {vertex.get_neighbours()};
    const auto& indirect_neighbours {vertex.get_indirect_neighbours()};

    do {
        // std::uniform_int_distribution<std::size_t> int_distribution(1, model_ref.model_params.vertices);
        std::uniform_int_distribution<std::size_t> int_distribution(1, model_ref.model_params.max_degree + 1);
        std::size_t drawn_colour {int_distribution(model_ref.rng)};

        auto match_colour {
            [=](const Vertex* const vertex) {
                return vertex->get_colour() == drawn_colour;
            }};

        auto it1 {std::find_if(
            neighbours.begin(),
            neighbours.end(),
            match_colour
        )};

        auto it2 {std::find_if(
            indirect_neighbours.begin(),
            indirect_neighbours.end(),
            match_colour
        )};


        if (it1 == neighbours.end() && it2 == indirect_neighbours.end()) {
            vertex.update_colour(drawn_colour);
            correctly_coloured = true;
        }
    } while (!correctly_coloured);
}

std::size_t Solver::find_first_avail_colour(Vertex& vertex) const
{
    const auto& neighbours {vertex.get_neighbours()};
    const auto& indirect_neighbours {vertex.get_indirect_neighbours()};

    std::set<std::size_t> forbidden_colours_set;
    for (const auto* const neighbour : neighbours)
        forbidden_colours_set.insert(neighbour->get_colour());

    for (const auto* const neighbour : indirect_neighbours)
        forbidden_colours_set.insert(neighbour->get_colour());

    std::vector<std::size_t> forbidden_colours(forbidden_colours_set.begin(), forbidden_colours_set.end());

    std::cout << "\n========================================";
    std::cout << "\n--> Vertex: " << vertex;
    std::cout << "\n----> Neighbours: [" << neighbours << "]";
    std::cout << "\n----> Indirect neighbours: [" << indirect_neighbours << "]";
    std::cout << "\n----> Forbidden set: [" << forbidden_colours_set << "]";
    std::cout << "\n----> Forbidden vector: [" << forbidden_colours << "]";
    std::cout << "\n========================================\n";

    std::size_t avail_colour {1};
    if (forbidden_colours.size() == 0) {
        return avail_colour;
    }

    if (forbidden_colours.front() > avail_colour)
        return avail_colour;

    for (std::size_t i {0}; i < forbidden_colours.size(); i++) {
        if (forbidden_colours.size() == 1) {
            return forbidden_colours.at(0) == 1 ? 2 : 1;
        }

        if (i == forbidden_colours.size() - 1) {
            return forbidden_colours.at(i) + 1;
        }

        if (forbidden_colours.at(i + 1) - forbidden_colours.at(i) > 1) {
            avail_colour = (forbidden_colours.at(i) + 1);
            break;
        }
    }

    return avail_colour;
}
