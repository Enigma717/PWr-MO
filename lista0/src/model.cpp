#include "../include/model.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <limits>

namespace
{
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

void Model::solve_knapsack()
{
    const std::size_t items_count {items.size()};
    const std::size_t capacity {model_params.capacity};
    std::vector<std::vector<int>> dp_matrix(items_count + 1, std::vector<int>(capacity + 1));

    for (std::size_t i = 0; i <= items_count; i++) {
        for (std::size_t w = 0; w <= capacity; w++) {
            if (i == 0 || w == 0)
                dp_matrix.at(i).at(w) = 0;
            else if (items.at(i - 1).weight <= w)
                dp_matrix.at(i).at(w) = std::max(items.at(i - 1).profit +
                    dp_matrix.at(i - 1).at(w - items.at(i - 1).weight), dp_matrix.at(i - 1).at(w));
            else
                dp_matrix.at(i).at(w) = dp_matrix.at(i - 1).at(w);
        }
    }

    int result {dp_matrix.at(items_count).at(capacity)};

    knapsack_value = result;
    knapsack_solution.resize(items_count);

    int weight = capacity;
    for (int i = items_count; i > 0 && result > 0; i--) {
        if (result != dp_matrix.at(i - 1).at(weight)) {
            knapsack_solution.at(items.at(i - 1).index - 1) = true;

            result = result - items.at(i - 1).profit;
            weight = weight - items.at(i - 1).weight;
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
        << "\n|-> Speed_to_weight_ratio: " << model_params.speed_to_weight_ratio
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

std::string Model::print_items() const
{
    std::stringstream stream;

    stream << "[";

    for (std::size_t i {0uz}; i < items.size(); i++) {
        const auto& item {items.at(i)};
        stream << "\nItem(" << item.index << "): [node_idx: "
            << item.node_index << ", profit: "
            << item.profit << ", weight: "
            << item.weight << ", ratio: "
            << item.ratio << "]";

        if (i != items.size() - 1)
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
        const std::size_t source {solution.at(i).index};
        const std::size_t destination {solution.at((i + 1) % solution.size()).index};
        const double distance {weights.at(source - 1).at(destination - 1)};

        objective_sum += distance;
    }

    return objective_sum;
}

double Model::evaluate_member_fitness(const Member& member)
{
    const std::vector<Node>& solution {member.solution};
    double objective_sum {static_cast<double>(knapsack_value)};

    for (std::size_t i {0uz}; i < solution.size(); i++) {
        const std::size_t source {solution.at(i).index};
        const std::size_t destination {solution.at((i + 1) % (solution.size())).index};
        const double distance {weights.at(source - 1).at(destination - 1)};

        // std::cout << "\nSource: " << source;
        // std::cout << "\nSource - 1: " << (source - 1);
        int item_weight {0};
        if (source != 1 && knapsack_solution.at(source - 2))
            item_weight = items.at(source - 2).weight;

        const double travel_time {calculate_travel_time(distance, item_weight)};

        objective_sum -= travel_time;
    }

    current_knapsack_weight = 0;

    return objective_sum;
}

double Model::calculate_travel_time(const double distance, const int item_weight)
{
    // std::cout << "\nKnapsack weight b4: " << knapsack_weight;
    current_knapsack_weight += item_weight;

    const double velocity {
        model_params.max_speed - (current_knapsack_weight * model_params.speed_to_weight_ratio)};

    // std::cout << "\nKnapsack weight after: " << knapsack_weight;
    // std::cout << "\nDistance: " << distance;
    // std::cout << "\nVelocity: " << velocity;
    // std::cout << "\nTTP distance: " << (distance / velocity) << "\n";

    return (distance / velocity);
}

std::vector<Node> Model::k_random_solution(std::size_t k_factor)
{
    std::vector<Node> solution {nodes};
    double distance {objective_function(solution)};
    double temp_distance {distance};

    for (std::size_t i {0uz}; i < k_factor; i++) {
        std::shuffle(solution.begin(), solution.end(), rng);
        temp_distance = objective_function(solution);

        if (distance > temp_distance)
            distance = temp_distance;
    }

    return solution;
}

std::vector<Node> Model::nearest_neighbour(std::size_t starting_node_index) const
{
    std::vector<bool> node_visit_status(model_params.dimension);
    std::vector<Node> solution {nodes};

    const std::size_t starting_node_position {starting_node_index - 1};
    const Node starting_node {nodes.at(starting_node_position)};

    solution.at(0) = starting_node;
    node_visit_status.at(starting_node_position) = true;
    std::size_t best_neighbour_position {0uz};

    for (std::size_t processed_position {0uz}; processed_position < model_params.dimension - 1;) {
        double best_distance {std::numeric_limits<double>::max()};

        for (std::size_t neighbour_position {0uz}; neighbour_position < node_visit_status.size(); neighbour_position++) {
            if (node_visit_status.at(neighbour_position) == false) {
                const std::size_t processed_node_position {solution.at(processed_position).index - 1};
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
