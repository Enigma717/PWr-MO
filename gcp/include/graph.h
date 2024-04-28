#pragma once

#include "./vertex.h"

#include <cstdint>
#include <vector>

class Graph {
public:
    Graph() = delete;
    Graph(std::size_t);

    void add_edge(std::size_t first_vertex_id, std::size_t second_vertex_id);

// private:
    std::vector<Vertex> vertices;
};
