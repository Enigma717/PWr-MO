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
    std::vector<bool> visited(vertices.size());

    bool bipartite_flag {bipartite_visit(&*vertices.begin(), visited)};

    if (bipartite_flag) {
        std::cout << "Graph is bipartite";

        std::vector<Vertex> first_set;
        std::vector<Vertex> second_set;

        auto find_color_0 = [](const Vertex& vertex)
        {
            return vertex.get_colour() == 0;
        };

        std::copy_if(
            vertices.begin(),
            vertices.end(),
            std::back_inserter(first_set),
            find_color_0);

        auto find_color_1 = [](const Vertex& vertex)
        {
            return vertex.get_colour() == 1;
        };

        std::copy_if(
            vertices.begin(),
            vertices.end(),
            std::back_inserter(second_set),
            find_color_1);

        std::cout << "\n |--> Vertices in first set: [" << first_set << "]";
        std::cout << "\n |--> Vertices in second set: [" << second_set << "]\n";

        return true;
    }
    else {
        std::cout << "Graph is not bipartite\n";

        return false;
    }
}

bool Graph::bipartite_visit(Vertex* vertex, std::vector<bool>& visited)
{
    if (const std::size_t id {vertex->get_id()}; id == 0) {
        visited.at(id) = true;
        vertex->update_colour(1);
    }

    for (const auto& neighbour : vertex->get_neighbours()) {
        const std::size_t neighbour_id {neighbour->get_id()};

        if (!visited.at(neighbour_id)) {
            visited.at(neighbour_id) = true;
            neighbour->update_colour(!vertex->get_colour());

            if (!bipartite_visit(neighbour, visited))
                return false;
        }
        else if (neighbour->get_colour() == vertex->get_colour()) {
            return false;
        }
    }

    return true;
}
