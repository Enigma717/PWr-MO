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
    void create_clusters();

    std::size_t variables_count {0uz};
    std::size_t matrix_size {0uz};

    std::vector<std::vector<double>> DSM;
    std::vector<std::vector<double>> DSM_entropy;
    std::vector<std::vector<double>> DSM_dist;
    std::vector<std::vector<std::size_t>> clusters;
    std::vector<std::vector<std::size_t>> unmerged_clusters;

private:
    void prepare_DSM();
    void add_cluster_to_matrix(
        const std::size_t first_index,
        const std::size_t second_index,
        const std::vector<std::size_t>& first_cluster,
        const std::vector<std::size_t>& second_cluster);
};
