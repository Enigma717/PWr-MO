#pragma once

#include "./node.h"

#include <cstdint>
#include <vector>

struct Member {
    std::vector<Node> solution;
    double fitness {0.0};
};
