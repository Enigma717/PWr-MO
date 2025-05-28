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

    void add_edge(const std::size_t source_id, const std::size_t destination_id);
    void reset_colouring();

    void BFS(const std::size_t starting_vertex_id);
    void DFS(const std::size_t starting_vertex_id);
    bool is_bipartite();

private:
    bool bipartite_visit(Vertex& vertex);

public:
    std::vector<Vertex> vertices;
    std::vector<std::size_t*> colours;
};
