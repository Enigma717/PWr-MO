#include "../include/model.h"
#include "../include/stream_operators.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>

namespace
{
    constexpr double optimum {7542.0};

    double euc_distance(Node first_node, Node second_node)
    {
        const double result {
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
    for (std::size_t i {0uz}; i < model_params.dimension; i++) {
        const Node& node {nodes.at(i)};

        for (std::size_t j {0uz}; j < model_params.dimension; j++) {
            const Node& neighbour {nodes.at(j)};

            if (node.index != neighbour.index)
                weights.at(i).at(j) = euc_distance(neighbour, node );
        }
    }
}

std::string Model::print_model_parms() const
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

std::string Model::print_nodes() const
{
    std::stringstream stream;

    stream << "[";

    for (std::size_t i {0uz}; i < nodes.size(); i++) {
        stream << nodes.at(i).index;

        if (i != nodes.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream.str();
}

std::string Model::print_weight_matrix() const
{
    std::stringstream stream;

    for (std::size_t i {0uz}; i < weights.size(); i++) {
        stream << "\nNode(" << i << ")\t[";

        for (std::size_t j {0uz}; j < weights.at(i).size(); j++) {
            const double neighbour {weights.at(i).at(j)};

            stream << neighbour;

            if (j != weights.at(i).size() - 1)
                stream << ", ";
        }

        stream << "]";
    }

    return stream.str();
}

double Model::objective_function(const std::vector<Node>& solution) const
{
    double objective_sum {0.0};

    for (std::size_t i {0uz}; i < solution.size(); i++) {
        const std::uint16_t source {solution.at(i).index};
        std::uint16_t destination {0};

        if (i == solution.size() - 1)
            destination = solution.at(0).index;
        else
            destination = solution.at(i + 1).index;

        const double distance {weights.at(source - 1).at(destination - 1)};

        objective_sum += distance;
    }

    return objective_sum;
}

double Model::prd(double objective_sum) const
{
    return 100.0 * ((objective_sum - optimum) / optimum);
}

std::vector<Node> Model::k_random_solution(std::uint64_t k)
{
    std::vector<Node> solution {nodes};
    double distance {objective_function(solution)};
    double temp_distance {distance};

    for (std::size_t i {0uz}; i < k; i++) {
        std::shuffle(solution.begin(), solution.end(), rng);
        temp_distance = objective_function(solution);

        if (distance > temp_distance)
            distance = temp_distance;
    }

    return solution;
}

std::vector<Node> Model::nearest_neighbour(std::uint16_t starting_node_index) const
{
    std::vector<bool> node_visit_status(model_params.dimension);
    std::vector<Node> solution {nodes};

    const std::uint16_t starting_node_position {
        static_cast<std::uint16_t>(starting_node_index - 1)};
    const Node starting_node {nodes.at(starting_node_position)};

    solution.at(0) = starting_node;
    node_visit_status.at(starting_node_position) = true;
    std::uint16_t best_neighbour_position {0u};

    for (std::size_t processed_position {0uz}; processed_position < model_params.dimension - 1;) {
        double best_distance {9999.0};

        for (std::size_t neighbour_position {0uz}; neighbour_position < node_visit_status.size(); neighbour_position++) {
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

        const Node neighbour_node {nodes.at(best_neighbour_position)};

        processed_position++;
        solution.at(processed_position) = neighbour_node;
        node_visit_status.at(best_neighbour_position) = true;
    }

    return solution;
}

std::vector<Node> Model::extended_nearest_neighbour() const
{
    std::vector<Node> solution {nodes};

    for (std::size_t starting_index {1uz}; starting_index <= model_params.dimension; starting_index++) {
        const std::vector<Node> temp_solution {nearest_neighbour(starting_index)};

        if (objective_function(solution) > objective_function(temp_solution))
            solution = temp_solution;
    }

    return solution;
}
