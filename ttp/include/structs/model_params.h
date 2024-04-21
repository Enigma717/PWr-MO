#pragma once

#include "../enums/problem_type.h"

#include <string>

struct ModelParams {
    std::string instance_name {""};
    ProblemType problem_type {ProblemType::bounded};
    std::size_t dimension {0};
    std::size_t num_of_items {0};
    std::size_t capacity {0};
    double min_speed {0.0};
    double max_speed {0.0};
    double renting_ratio {0.0};
    double speed_to_weight_ratio {0.0};
};
