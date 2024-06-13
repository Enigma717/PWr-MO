#pragma once

#include <string>

struct ModelParams {
    std::string instance_name {""};
    std::size_t vertices {0};
    std::size_t edges {0};
    std::size_t max_degree {0};
    std::size_t optimum {0};
};
