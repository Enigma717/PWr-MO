#include "linkage_tree_builder.hpp"
#include "utility_operators.hpp"

#include <iostream>
#include <map>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

LinkageTreeBuilder::LinkageTreeBuilder(std::size_t variables_count)
: variables_count {variables_count},
  matrix_size {variables_count}
{}

void LinkageTreeBuilder::print_DSM() const
{
    std::cout << "\nDSM:";
    for (const auto& row : DSM) {
        std::cout << "\n|";
        for (const auto value : row) {
            std::cout << "  " << std::setw(7) << std::fixed << std::setprecision(4) << value;
        }
        std::cout << "  |";
    }
}

void LinkageTreeBuilder::print_DSM_entropy() const
{
    std::cout << "\nDSM entropy:";
    for (const auto& row : DSM_entropy) {
        std::cout << "\n|";
        for (const auto value : row) {
            std::cout << "  " << std::setw(7) << std::fixed << std::setprecision(4) << value;
        }
        std::cout << "  |";
    }
}

void LinkageTreeBuilder::print_DSM_dist() const
{
    std::cout << "\nDSM distance:";
    for (const auto& row : DSM_dist) {
        std::cout << "\n|";
        for (const auto value : row) {
            std::cout << "  " << std::setw(7) << std::fixed << std::setprecision(4) << value;
        }
        std::cout << "  |";
    }
}

void LinkageTreeBuilder::prepare_DSM()
{
    matrix_size = variables_count;

    DSM.clear();
    DSM.resize(matrix_size);

    for (auto& row : DSM) {
        row.clear();
        row.resize(matrix_size);
    }

    DSM_entropy.clear();
    DSM_entropy.resize(matrix_size);

    for (auto& row : DSM_entropy) {
        row.clear();
        row.resize(matrix_size);
    }

    DSM_dist.clear();
    DSM_dist.resize(matrix_size);

    for (auto& row : DSM_dist) {
        row.clear();
        row.resize(matrix_size);
    }
}

void LinkageTreeBuilder::calculate_DSM(const std::vector<Solution>& subpopulation)
{
    prepare_DSM();

    const auto subpopulation_size {subpopulation.size()};

    for (auto x {0uz}; x < matrix_size; x++) {
        for (auto y {x}; y < matrix_size; y++) {
            if (x == y) {
                DSM[x][y] = -1.0;
                DSM_entropy[x][y] = -1.0;
                DSM_dist[x][y] = 9.0;
                continue;
            }

            std::map<std::size_t, std::size_t> first_col_colour_counts;
            std::map<std::size_t, std::size_t> second_col_colour_counts;
            std::map<ColoursPair, std::size_t> colors_pairs_counts;

            std::map<std::size_t, double> first_column_colour_probs;
            std::map<std::size_t, double> second_column_colour_probs;
            std::map<ColoursPair, double> colors_pairs_probs;

            for (auto i {0uz}; i < subpopulation_size; i++) {
                const auto& vertices {subpopulation.at(i).graph.vertices};
                const auto first_colour {vertices.at(x).colour};
                const auto second_colour {vertices.at(y).colour};
                const auto colours_pair {std::make_tuple(first_colour, second_colour)};

                auto first_it(first_col_colour_counts.find(first_colour));
                if (first_it != first_col_colour_counts.end())
                    first_it->second++;
                else
                    first_col_colour_counts[first_colour] = 1;

                auto second_it(second_col_colour_counts.find(second_colour));
                if (second_it != second_col_colour_counts.end())
                    second_it->second++;
                else
                second_col_colour_counts[second_colour] = 1;

                auto pair_it(colors_pairs_counts.find(colours_pair));
                if (pair_it != colors_pairs_counts.end())
                    pair_it->second++;
                else
                    colors_pairs_counts[colours_pair] = 1;
            }

            for (const auto& elem : first_col_colour_counts) {
                const double probability {
                    static_cast<double>(elem.second) /
                    static_cast<double>(subpopulation_size)};
                first_column_colour_probs[elem.first] = probability;
            }

            for (const auto& elem : second_col_colour_counts) {
                const double probability {
                    static_cast<double>(elem.second) /
                    static_cast<double>(subpopulation_size)};
                second_column_colour_probs[elem.first] = probability;
            }

            for (const auto& elem : colors_pairs_counts) {
                const double probability {
                    static_cast<double>(elem.second) /
                    static_cast<double>(subpopulation_size)};
                colors_pairs_probs[elem.first] = probability;
            }

            for (const auto& elems : colors_pairs_counts) {
                const auto first_color {elems.first.first};
                const auto second_color {elems.first.second};
                const auto pair_prob {colors_pairs_probs[elems.first]};
                const auto first_column_prob {first_column_colour_probs[first_color]};
                const auto second_column_prob {second_column_colour_probs[second_color]};

                const double MI {pair_prob * std::log(pair_prob / (first_column_prob * second_column_prob))};
                const double entropy {-pair_prob * std::log(pair_prob)};

                DSM[x][y] += MI;
                DSM[y][x] += MI;

                DSM_entropy[x][y] += entropy;
                DSM_entropy[y][x] += entropy;
            }

            const double distance_or_nan {(DSM_entropy[x][y] - DSM[x][y]) / DSM_entropy[x][y]};
            const double distance {std::isnan(distance_or_nan) ? 0.0 : distance_or_nan};

            DSM_dist[x][y] = distance;
            DSM_dist[y][x] = distance;
        }
    }
}

void LinkageTreeBuilder::create_clusters()
{
    clusters.clear();
    unmerged_clusters.clear();

    for (std::size_t i {0uz}; i < variables_count; i++) {
        clusters.push_back({i});
        unmerged_clusters.push_back({i});
    }

    while (unmerged_clusters.size() > 1) {
        double smallest_value {DSM_dist[0][1]};
        std::vector<std::size_t> smallest_value_indexes {0, 1};
        for (std::size_t i {0uz}; i < matrix_size - 1; i++) {
            const auto curr_smallest_value_it {std::min_element(DSM_dist[i].begin() + (i + 1), DSM_dist[i].end())};

            if (*curr_smallest_value_it < smallest_value) {
                smallest_value = *curr_smallest_value_it;
                smallest_value_indexes[0] = i;
                smallest_value_indexes[1] = curr_smallest_value_it - DSM_dist[i].begin();
            }
        }

        const auto first_index {smallest_value_indexes[0]};
        const auto second_index {smallest_value_indexes[1]};
        const auto first_cluster {unmerged_clusters[first_index]};
        const auto second_cluster {unmerged_clusters[second_index]};


        std::vector<std::size_t> new_cluster;
        new_cluster.reserve(first_cluster.size() + second_cluster.size() ); // preallocate memory
        new_cluster.insert(new_cluster.end(), first_cluster.begin(), first_cluster.end());
        new_cluster.insert(new_cluster.end(), second_cluster.begin(), second_cluster.end());
        std::sort(new_cluster.begin(), new_cluster.end());

        unmerged_clusters.push_back(new_cluster);
        std::erase(unmerged_clusters, std::vector<std::size_t>{first_cluster.begin(), first_cluster.end()});
        std::erase(unmerged_clusters, std::vector<std::size_t>{second_cluster.begin(), second_cluster.end()});
        clusters.push_back(new_cluster);

        if (smallest_value <= 1e-7f) {
            std::erase(clusters, std::vector<std::size_t>{first_cluster.begin(), first_cluster.end()});
            std::erase(clusters, std::vector<std::size_t>{second_cluster.begin(), second_cluster.end()});
        }

        add_cluster_to_matrix(first_index, second_index, first_cluster, second_cluster);

        DSM_dist.erase(DSM_dist.begin() + second_index);
        for (auto& row : DSM_dist)
            row.erase(row.begin() + second_index);

        DSM_dist.erase(DSM_dist.begin() + first_index);
        for (auto& row : DSM_dist)
            row.erase(row.begin() + first_index);
    }

    clusters.pop_back();
    std::sort(
        clusters.begin(),
        clusters.end(),
        [&](const auto& cluster1, const auto& cluster2){
            return cluster1.size() < cluster2.size();
        });
}

void LinkageTreeBuilder::add_cluster_to_matrix(
    const std::size_t first_index,
    const std::size_t second_index,
    const std::vector<std::size_t>& first_cluster,
    const std::vector<std::size_t>& second_cluster)
{
    std::vector<double> new_distances;
    new_distances.reserve(matrix_size + 1);

    const auto quotient {first_cluster.size() + second_cluster.size()};
    const auto first_factor {static_cast<double>(first_cluster.size()) / quotient};
    const auto second_factor {static_cast<double>(second_cluster.size()) / quotient};

    for (std::size_t i {0uz}; i < matrix_size; i++) {
        auto& row {DSM_dist[i]};
        const auto first_dist {row[first_index]};
        const auto second_dist {row[second_index]};
        const auto new_dist {(first_factor * first_dist) + (second_factor * second_dist)};

        row.push_back(new_dist);
        new_distances.push_back(new_dist);
    }

    new_distances.push_back(9.0);
    DSM_dist.push_back(new_distances);

    matrix_size--;
}
