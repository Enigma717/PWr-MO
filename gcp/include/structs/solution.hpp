#pragma once

#include "graph.hpp"

struct Solution {
    Solution() = default;
    Solution(Graph&& graph) : graph {std::move(graph)} {};

    bool operator==(const Solution&) const = default;

    Graph graph;
    std::size_t fitness {0};
};

template<>
struct std::hash<Solution>
{
    std::size_t operator()(const Solution& solution) const noexcept
    {
        const auto& colours {solution.graph.colours};
        std::size_t seed {colours.size()};

        for(const auto* colour : colours)
            seed ^= *colour + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};
