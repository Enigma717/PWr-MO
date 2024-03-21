#include "../include/model.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

namespace
{
    constexpr double optimum {7542.0};

    double euc_distance(Node first_node, Node second_node)
    {
        const auto result {
            std::sqrt(std::pow((first_node.x_coord - second_node.x_coord), 2) +
            std::pow((first_node.y_coord - second_node.y_coord), 2))};

        return result;
    }
}


Model::Model() : loader{*this}, genetic_solver{*this} {}

void Model::load_file(const std::string& file_path)
{
    loader.parse_instance(file_path);
}

void Model::create_weight_matrix()
{
    for (auto i {0}; i < model_params.dimension; i++) {
        const auto& node {nodes.at(i)};

        for (auto j {0}; j < model_params.dimension; j++) {
            const auto& neighbour {nodes.at(j)};

            if (node.index != neighbour.index)
                weights.at(i).at(j) = euc_distance(neighbour, node );
        }
    }
}

std::string Model::print_model_parms()
{
    std::stringstream stream;
    stream << "-> Model parameters:"
        << "\n|-> Name: " << model_params.instance_name
        << "\n|-> Problem_type: " << model_params.problem_type
        << "\n|-> Dimension: " << model_params.dimension
        << "\n|-> Num_of_items: " << model_params.num_of_items
        << "\n|-> Capacity: " << model_params.capacity
        << "\n|-> Min_speed: " << model_params.min_speed
        << "\n|-> Max_speed: " << model_params.max_speed
        << "\n|-> Renting_ratio: " << model_params.renting_ratio
        << "\n|-> Nodes size: " << nodes.size()
        << "\n\\-> Items size: " << items.size();

    return stream.str();
}

std::string Model::print_nodes()
{
    std::stringstream stream;

    stream << "[";

    for (std::size_t i {0}; i < nodes.size(); i++) {
        // stream << i << ":" << nodes.at(i).index;
        stream << nodes.at(i).index;

        if (i != nodes.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream.str();
}

std::string Model::print_nodes(const std::vector<Node>& solution)
{
    std::stringstream stream;

    stream << "[";

    for (std::size_t i {0}; i < solution.size(); i++) {
        // stream << i << ":" << solution.at(i).index;
        stream << solution.at(i).index;

        if (i != solution.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream.str();
}

std::string Model::print_weight_matrix()
{
    std::stringstream stream;

    for (std::size_t i {0}; i < weights.size(); i++) {
        stream << "\nNode(" << i << ")\t[";

        for (std::size_t j {0}; j < weights.at(i).size(); j++) {
            const auto neighbour {weights.at(i).at(j)};

            // stream << j << ":" << neighbour;
            stream << neighbour;

            if (j != weights.at(i).size() - 1)
                stream << ", ";
        }

        stream << "]";
    }

    return stream.str();
}

double Model::objective_function()
{
    double objective_sum {0.0};

    for (std::size_t i {0}; i < nodes.size(); i++) {
        const auto source {nodes.at(i).index};
        auto destination {0};

        if (i == nodes.size() - 1)
            destination = nodes.at(0).index;
        else
            destination = nodes.at(i + 1).index;

        const auto distance {weights.at(source - 1).at(destination - 1)};

        objective_sum += distance;

        // std::cout << ">>> SOURCE: " << source << "\n";
        // std::cout << ">>> DESTINATION: " << destination << "\n";
        // std::cout << ">>> DIST: " << distance << "\n";
        // std::cout << ">>> SUM: " << objective_sum << "\n";
    }

    return objective_sum;
}

double Model::objective_function(const std::vector<Node>& solution)
{
    double objective_sum {0.0};

    for (std::size_t i {0}; i < solution.size(); i++) {
        const auto source {solution.at(i).index};
        auto destination {0};

        if (i == solution.size() - 1)
            destination = solution.at(0).index;
        else
            destination = solution.at(i + 1).index;

        const auto distance {weights.at(source - 1).at(destination - 1)};

        objective_sum += distance;
    }

    return objective_sum;
}

double Model::prd(double objective_sum)
{
    return 100.0 * ((objective_sum - optimum) / optimum);
}

std::vector<Node> Model::k_random_solution(std::uint64_t k)
{
    std::vector<Node> solution {nodes};
    double distance {objective_function(solution)};
    double temp_distance {distance};

    for (size_t i {0}; i < k; i++) {
        std::shuffle(solution.begin(), solution.end(), rng);
        temp_distance = objective_function(solution);

        if (distance > temp_distance)
            distance = temp_distance;
    }

    return solution;
}

std::vector<Node> Model::nearest_neighbour(std::uint16_t starting_node_index)
{
    std::vector<bool> node_visit_status(model_params.dimension);
    std::vector<Node> solution {nodes};

    const std::uint16_t starting_node_position {
        static_cast<std::uint16_t>(starting_node_index - 1)};
    const Node starting_node {nodes.at(starting_node_position)};

    solution.at(0) = starting_node;
    node_visit_status.at(starting_node_position) = true;
    std::uint16_t best_neighbour_position {0};

    for (size_t processed_position {0}; processed_position < model_params.dimension - 1;) {
        double best_distance {9999.0};

        for (size_t neighbour_position {0}; neighbour_position < node_visit_status.size(); neighbour_position++) {
            if (node_visit_status.at(neighbour_position) == false) {
                const std::uint16_t processed_node_position {
                    static_cast<std::uint16_t>(solution.at(processed_position).index - 1)};
                double temp_distance {weights.at(processed_node_position).at(neighbour_position)};

                if (best_distance > temp_distance) {
                    best_neighbour_position = neighbour_position;
                    best_distance = temp_distance;
                }
            }
        }

        const auto neighbour_node {nodes.at(best_neighbour_position)};

        processed_position++;
        solution.at(processed_position) = neighbour_node;
        node_visit_status.at(best_neighbour_position) = true;
    }

    return solution;
}

std::vector<Node> Model::extended_nearest_neighbour()
{
    std::vector<Node> solution {nodes};

    for (size_t starting_index {1}; starting_index <= model_params.dimension; starting_index++) {
        const auto temp_solution {nearest_neighbour(starting_index)};
        if (objective_function(solution) > objective_function(temp_solution))
            solution = temp_solution;
    }

    return solution;
}
