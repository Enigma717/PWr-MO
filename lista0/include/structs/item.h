#pragma once

#include <cstdint>

struct Item {
    std::size_t index {0uz};
    std::size_t node_index {0uz};
    int profit {0};
    int weight {0};
    double ratio {0.0};
};
