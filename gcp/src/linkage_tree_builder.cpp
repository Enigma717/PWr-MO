#include "linkage_tree_builder.hpp"

#include <iostream>
#include <map>
#include <iomanip>
#include <cmath>

LinkageTreeBuilder::LinkageTreeBuilder(std::size_t variables_count)
: variables_count {variables_count}
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
    DSM.clear();
    DSM.resize(variables_count);

    for (auto& row : DSM) {
        row.clear();
        row.resize(variables_count);
    }

    DSM_entropy.clear();
    DSM_entropy.resize(variables_count);

    for (auto& row : DSM_entropy) {
        row.clear();
        row.resize(variables_count);
    }

    DSM_dist.clear();
    DSM_dist.resize(variables_count);

    for (auto& row : DSM_dist) {
        row.clear();
        row.resize(variables_count);
    }
}

void LinkageTreeBuilder::calculate_DSM(const std::vector<Solution>& subpopulation)
{
    prepare_DSM();

    const auto subpopulation_size {subpopulation.size()};

    for (auto x {0uz}; x < variables_count; x++) {
        for (auto y {x}; y < variables_count; y++) {
            if (x == y) {
                DSM[x][y] = -1.0;
                DSM_entropy[x][y] = -1.0;
                DSM_dist[x][y] = -1.0;
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

                // std::cout << "\nIND" << i << "[" << x << "]: " << first_colour << " | IND" << i << "[" << y << "]: " << second_colour;

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

            // std::cout << "\nFirst column counts:";
            // for (const auto& elems : first_col_colour_counts) {
                // std::cout << "\n[" << elems.first << "]: " << elems.second;
            // }

            // std::cout << "\nSecond column counts:";
            // for (const auto& elems : second_col_colour_counts) {
                // std::cout << "\n[" << elems.first << "]: " << elems.second;
            // }

            // std::cout << "\nColors pairs counts:";
            // for (const auto& elems : colors_pairs_counts) {
            //     std::cout << "\n[" << elems.first.first << ", " << elems.first.second << "]: " << elems.second;
            // }

            // std::cout << "\nFirst column probabilities:";
            // for (const auto& elems : first_column_colour_probs) {
            //     std::cout << "\np_" << x << "[" << elems.first << "]: " << elems.second;
            // }

            // std::cout << "\nSecond column probabilities:";
            // for (const auto& elems : second_column_colour_probs) {
            //     std::cout << "\np_" << y << "[" << elems.first << "]: " << elems.second;
            // }

            // std::cout << "\nColors pairs probabilities:";
            // for (const auto& elems : colors_pairs_probs) {
            //     std::cout << "\np_" << x << y << "[" << elems.first.first << ", " << elems.first.second << "]: " << elems.second;
            // }

            for (const auto& elems : colors_pairs_counts) {
                const auto first_color {elems.first.first};
                const auto second_color {elems.first.second};
                const auto pair_prob {colors_pairs_probs[elems.first]};
                const auto first_column_prob {first_column_colour_probs[first_color]};
                const auto second_column_prob {second_column_colour_probs[second_color]};

                const double MI {pair_prob * std::log(pair_prob / (first_column_prob * second_column_prob))};
                const double entropy {-pair_prob * std::log(pair_prob)};

                // std::cout << "\nMI: p_" << x << y << "[" << first_color << ", " << second_color << "]"
                //     << " * ln(p_" << x << y << "[" << first_color << ", " << second_color << "]"
                //     << " / (p_" << x << "[" << first_color << "] * p_" << y << "[" << second_color << "]))"
                //     << " = " << pair_prob << " * ln(" << pair_prob << " / (" << first_column_prob << " * " << second_column_prob  << "))"
                //     << " = " << MI;
                // std::cout << "\nH: -p_" << x << y << "[" << first_color << ", " << second_color << "]"
                //     << " * ln(p_" << x << y << "[" << first_color << ", " << second_color << "])"
                //     << " / (p_" << x << "[" << first_color << "] * p_" << y << "[" << second_color << "]))"
                //     << " = " << -pair_prob << " * ln(" << pair_prob  << ")"
                //     << " = " << entropy;

                DSM[x][y] += MI;
                DSM[y][x] += MI;

                DSM_entropy[x][y] += entropy;
                DSM_entropy[y][x] += entropy;
            }

            const double distance {(DSM_entropy[x][y] - DSM[x][y]) / DSM_entropy[x][y]};

            // std::cout << "\nD: (H[" << x << ", " << y << "]"
            //     << " - MI[" << x << ", " << y << "])"
            //     << " / H[" << x << ", " << y << "]"
            //     << " = " << "(" << DSM_entropy[x][y] << " - " << DSM[x][y] << ") / " << DSM_entropy[x][y]
            //     << " = " << distance;

            DSM_dist[x][y] = distance;
            DSM_dist[y][x] = distance;
        }
    }

    // print_DSM();
    // print_DSM_entropy();
    // print_DSM_dist();
}