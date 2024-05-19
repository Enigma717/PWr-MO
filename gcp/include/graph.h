#pragma once

#include "./vertex.h"

#include <cstdint>
#include <vector>

class Graph {
public:
    Graph() = delete;
    Graph(const std::size_t);

    void add_edge(const std::size_t source_id, const std::size_t destination_id);

    void BFS(const std::size_t starting_vertex_id);
    void DFS(const std::size_t starting_vertex_id);
    bool is_bipartite();

private:
    bool bipartite_visit(Vertex* vertex, std::vector<bool>& visited);

    std::vector<Vertex> vertices;
    std::vector<std::size_t*> colours;
};
