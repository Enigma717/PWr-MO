#pragma once

#include "../enums/problem_type.h"

#include <string>

struct ModelParams {
    std::string instance_name {""};
    ProblemType problem_type {ProblemType::bounded};
    int dimension {0};
    int num_of_items {0};
    int capacity {0};
    double min_speed {0.0};
    double max_speed {0.0};
    double renting_ratio {0.0};
};
