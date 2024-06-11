#pragma once

#include "../graph.h"

struct Solution {
    Solution() = default;
    Solution(Graph&& graph) : graph{std::move(graph)} {};

    Graph graph;
    std::size_t fitness {0};
};
