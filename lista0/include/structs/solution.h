#pragma once

#include "./node.h"
#include "./item.h"

#include <vector>

struct Solution {
    std::vector<Node> route;
    std::vector<Item> penalized_items;
    std::vector<bool> packing_plan;
    std::size_t knapsack_weight {0uz};
    int knapsack_value {0};
    double fitness {0.0};
};
