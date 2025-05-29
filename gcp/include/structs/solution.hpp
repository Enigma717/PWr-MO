#pragma once

#include "graph.hpp"

struct Solution {
    Solution() = default;
    Solution(Graph&& graph) : graph {std::move(graph)} {};

    Graph graph;
    std::size_t fitness {0};
};
