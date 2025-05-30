#pragma once

#include "structs/solution.hpp"

#include <vector>
#include <map>

class LinkageTreeBuilder {
public:
    using ColoursPair = std::pair<std::size_t, std::size_t>;

    LinkageTreeBuilder(std::size_t variables_count);

    void print_DSM() const;
    void print_DSM_entropy() const;
    void print_DSM_dist() const;
    void calculate_DSM(const std::vector<Solution>& subpopulation);

    std::size_t variables_count {0uz};

    std::vector<std::vector<double>> DSM;
    std::vector<std::vector<double>> DSM_entropy;
    std::vector<std::vector<double>> DSM_dist;

private:
    void prepare_DSM();
};
