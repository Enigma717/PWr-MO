#include "../include/graph.h"
#include "../include/utility_operators.h"

#include <iostream>

Graph::Graph(std::size_t size)
{
    vertices.reserve(size);

    for (std::size_t i {0}; i < size; i++)
        vertices.emplace_back();

}

void Graph::add_edge(std::size_t first_vertex_id, std::size_t second_vertex_id)
{
    Vertex& first_vertex {vertices[first_vertex_id]};
    Vertex& second_vertex {vertices[second_vertex_id]};

    first_vertex.update_neighbourship(second_vertex);
    second_vertex.update_neighbourship(first_vertex);
}
