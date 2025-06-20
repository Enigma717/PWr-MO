#pragma once

#include "vertex.hpp"

#include <cstdint>
#include <vector>

class Graph {
public:
    Graph() = default;
    ~Graph() = default;
    Graph(const Graph&);
    Graph& operator=(const Graph&);
    Graph(Graph&&) = default;
    Graph& operator=(Graph&&) = default;
    Graph(const std::size_t);

    bool operator==(const Graph&) const = default;

    void add_edge(const std::size_t source_id, const std::size_t destination_id);
    void reset_colouring();

    void BFS(const std::size_t starting_vertex_id);
    void DFS(const std::size_t starting_vertex_id);
    bool is_bipartite();
    bool bipartite_visit(Vertex& vertex);

    std::vector<Vertex> vertices;
    std::vector<std::size_t*> colours;
};
