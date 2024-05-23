#include "../include/graph.h"
#include "../include/utility_operators.h"

#include <stack>
#include <queue>
#include <iostream>
#include <algorithm>
#include <iterator>

Graph::Graph(const std::size_t size)
{
    vertices.reserve(size);
    colours.reserve(size);

    for (std::size_t i {0}; i < size; i++) {
        vertices.emplace_back();
        colours.push_back(&(vertices.at(i).get_colour_ref()));
    }
}

void Graph::add_edge(const std::size_t source_id, const std::size_t destination_id)
{
    Vertex& source {vertices.at(source_id)};
    Vertex& destination {vertices.at(destination_id)};

    source.update_neighbourship_with(destination);
}

void Graph::reset_colouring()
{
    for (auto& vertex : vertices) {
        vertex.update_colour(0);
    }
}

void Graph::BFS(const std::size_t starting_vertex_id)
{
    std::vector<bool> visited(vertices.size());
    std::queue<Vertex*> queue;

    visited.at(starting_vertex_id) = true;
    queue.push(&vertices.at(starting_vertex_id));

    while (!queue.empty()) {
        const Vertex* current {queue.front()};
        queue.pop();

        std::cout << *current << " -> ";

        for (const auto& neighbour : current->get_neighbours()) {
            const std::size_t neighbour_id {neighbour->get_id()};

            if (!visited.at(neighbour_id)) {
                visited.at(neighbour_id) = true;
                queue.push(neighbour);
            }
        }
    }

    std::cout << "|\n";
}

void Graph::DFS(const std::size_t starting_vertex_id)
{
    std::vector<bool> visited(vertices.size());
    std::stack<Vertex*> stack;

    visited.at(starting_vertex_id) = true;
    stack.push(&vertices.at(starting_vertex_id));

    while (!stack.empty()) {
        const Vertex* current {stack.top()};
        stack.pop();

        std::cout << *current << " -> ";

        for (const auto& neighbour : current->get_neighbours()) {
            const std::size_t neighbour_id {neighbour->get_id()};

            if (!visited.at(neighbour_id)) {
                visited.at(neighbour_id) = true;
                stack.push(neighbour);
            }
        }
    }

    std::cout << "|\n";
}

bool Graph::is_bipartite()
{
    for (auto& vertex : vertices) {
        if (vertex.get_colour() == 0) {
            if (!bipartite_visit(vertex)) {
                std::cout << "Graph is not bipartite\n";

                return false;
            }
        }
    }

    std::cout << "Graph is bipartite";

    std::vector<Vertex> first_set;
    std::vector<Vertex> second_set;

    auto find_first_colour = [](const Vertex& vertex)
    {
        return vertex.get_colour() == 1;
    };

    std::copy_if(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(first_set),
        find_first_colour);

    auto find_second_colour = [](const Vertex& vertex)
    {
        return vertex.get_colour() == 2;
    };

    std::copy_if(
        vertices.begin(),
        vertices.end(),
        std::back_inserter(second_set),
        find_second_colour);

    std::cout << "\n |--> Vertices in first set: [" << first_set << "]";
    std::cout << "\n |--> Vertices in second set: [" << second_set << "]\n";

    return true;
}

bool Graph::bipartite_visit(Vertex& vertex)
{
    std::queue<Vertex*> queue;
    vertex.update_colour(1);
    queue.push(&vertex);

    while (!queue.empty()) {
        const Vertex* const current {queue.front()};
        queue.pop();

        for (const auto& neighbour : current->get_neighbours()) {
            if (neighbour->get_colour() == 0) {
                neighbour->update_colour(3 - current->get_colour());
                queue.push(neighbour);
            }
            else if (neighbour->get_colour() == current->get_colour())
                return false;
        }
    }

    return true;
}
