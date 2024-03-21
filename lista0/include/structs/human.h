#pragma once

#include "./node.h"

#include <cstdint>
#include <vector>

struct Human {
    std::vector<Node> solution;
    double objective {0.0};
};
