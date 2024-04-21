#include "../include/model.h"
#include "../include/utility_operators.h"
#include "../include/structs/solution.h"

#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <iomanip>

namespace
{
    constexpr double maximum_load_threshold {0.67};
    constexpr double initial_sa_temperature {1'000'000.0};
    constexpr double cooling_rate_value {0.9999};
    constexpr std::size_t maximum_sa_iterations {1'000'000uz};

    double euc_distance(Node first_node, Node second_node)
    {
        const double result {
            std::sqrt(std::pow((first_node.x_coord - second_node.x_coord), 2) +
            std::pow((first_node.y_coord - second_node.y_coord), 2))};

        return result;
    }

    constexpr double ln_factor(const std::size_t iteration, const std::size_t dimension)
    {
        return (std::log(iteration) / std::log(dimension));
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
        << "\n|-> Speed_to_weight_ratio: " << model_params.speed_to_weight_ratio
        << "\n|-> Nodes size: " << nodes.size()
        << "\n\\-> Items size: " << items.size();

    return stream.str();
}

std::string Model::print_nodes() const
{
    std::stringstream stream;

    stream << nodes;

    return stream.str();
}

std::string Model::print_items() const
{
    std::stringstream stream;

    stream << items;

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

double Model::objective_function(const std::vector<Node>& route) const
{
    double objective_sum {0.0};

    for (std::size_t i {0uz}; i < route.size(); i++) {
        const std::size_t source {route.at(i).index};
        const std::size_t destination {route.at((i + 1) % route.size()).index};
        const double distance {weights.at(source - 1).at(destination - 1)};

        objective_sum += distance;
    }

    return objective_sum;
}

double Model::evaluate_solution_fitness(Solution& solution)
{
    const std::vector<Node>& route {solution.route};
    const std::vector<bool>& packing_plan {solution.packing_plan};
    std::size_t& knapsack_weight {solution.knapsack_weight};
    double objective_sum {static_cast<double>(solution.knapsack_value)};

    for (std::size_t i {0uz}; i < route.size(); i++) {
        const std::size_t source {route.at(i).index};
        const std::size_t destination {route.at((i + 1) % (route.size())).index};
        const double distance {weights.at(source - 1).at(destination - 1)};

        int item_weight {0};
        if (source != 1 && packing_plan.at(source - 2))
            item_weight = items.at(source - 2).weight;

        const double travel_time {calculate_travel_time(
            distance, item_weight, knapsack_weight)};

        objective_sum -= travel_time;
    }

    return objective_sum;
}

double Model::calculate_travel_time(const double distance, const int item_weight, std::size_t& knapsack_weight)
{
    knapsack_weight += item_weight;

    const double velocity {
        model_params.max_speed - (knapsack_weight * model_params.speed_to_weight_ratio)};

    return (distance / velocity);
}

Solution Model::k_random_solution(std::size_t k_factor)
{
    std::vector<Node> route {nodes};
    Solution best_solution;

    best_solution.fitness = -999'999'999.0;

    for (std::size_t i {0uz}; i < k_factor; i++) {
        std::shuffle(route.begin(), route.end(), rng);

        Solution temp_solution;
        temp_solution.route = route;
        temp_solution.penalized_items = penalize_item_values(route);
        temp_solution.packing_plan = solve_knapsack_greedy(temp_solution);
        temp_solution.fitness = evaluate_solution_fitness(temp_solution);

        if (best_solution < temp_solution)
            best_solution = temp_solution;
    }

    return best_solution;
}

std::vector<Node> Model::nearest_neighbour(std::size_t starting_node_index)
{
    std::vector<bool> node_visit_status(model_params.dimension);
    std::vector<Node> route {nodes};

    const std::size_t starting_node_position {starting_node_index - 1};
    const Node starting_node {nodes.at(starting_node_position)};

    route.at(0) = starting_node;
    node_visit_status.at(starting_node_position) = true;
    std::size_t best_neighbour_position {0uz};

    for (std::size_t processed_position {0uz}; processed_position < model_params.dimension - 1;) {
        double best_distance {std::numeric_limits<double>::max()};

        for (std::size_t neighbour_position {0uz}; neighbour_position < node_visit_status.size(); neighbour_position++) {
            if (node_visit_status.at(neighbour_position) == false) {
                const std::size_t processed_node_position {route.at(processed_position).index - 1};
                double temp_distance {weights.at(processed_node_position).at(neighbour_position)};

                if (best_distance > temp_distance) {
                    best_neighbour_position = neighbour_position;
                    best_distance = temp_distance;
                }
            }
        }

        const Node neighbour_node {nodes.at(best_neighbour_position)};

        processed_position++;
        route.at(processed_position) = neighbour_node;
        node_visit_status.at(best_neighbour_position) = true;
    }

    return route;
}

Solution Model::extended_nearest_neighbour()
{
    std::vector<Node> route {nodes};
    Solution best_solution;

    best_solution.fitness = -999'999'999.0;

    for (std::size_t starting_index {1uz}; starting_index <= model_params.dimension; starting_index++) {
        const std::vector<Node> temp_route {nearest_neighbour(starting_index)};

        Solution temp_solution;
        temp_solution.route = temp_route;
        temp_solution.penalized_items = penalize_item_values(temp_route);
        temp_solution.packing_plan = solve_knapsack_greedy(temp_solution);
        temp_solution.fitness = evaluate_solution_fitness(temp_solution);

        if (best_solution < temp_solution)
            best_solution = temp_solution;
    }

    return best_solution;
}

Solution Model::simulated_annealing()
{
    const double cooling_rate {cooling_rate_value};
    double current_temperature {initial_sa_temperature};

    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    Solution best_solution;
    best_solution.route = nodes;
    best_solution.penalized_items = penalize_item_values(best_solution.route);
    best_solution.packing_plan = solve_knapsack_greedy(best_solution);
    best_solution.fitness = evaluate_solution_fitness(best_solution);

    for (std::size_t iteration {0}; iteration < maximum_sa_iterations; iteration++) {
        Solution neighbour;
        neighbour.route = process_invert_mutation(best_solution.route);
        neighbour.penalized_items = penalize_item_values(neighbour.route);
        neighbour.packing_plan = solve_knapsack_greedy(neighbour);
        neighbour.fitness = evaluate_solution_fitness(neighbour);

        const double fitness_diff {neighbour.fitness - best_solution.fitness};
        const double probability {distribution(rng)};
        const double exp {std::exp(fitness_diff / current_temperature)};

        if (fitness_diff >= 0) {
            best_solution = neighbour;
        }
        else {
            if (probability < exp)
                best_solution = neighbour;
        }

        current_temperature *= cooling_rate;
    }

    return best_solution;
}

Solution Model::genetic_algorithm()
{
    return genetic_solver.solve();
}


PackingPlan Model::solve_knapsack_greedy(Solution& solution)
{
    const std::size_t items_count {items.size()};
    const std::size_t capacity {model_params.capacity};

    std::sort(
        solution.penalized_items.begin(),
        solution.penalized_items.end(),
        [](Item& item1, Item& item2) { return item1.ratio > item2.ratio; });

    std::vector<bool> knapsack_solution(items_count);
    std::size_t current_weight {0uz};
    int current_profit {0};

    for (std::size_t i {0}; i < items_count; i++) {
        const auto& current_item {solution.penalized_items.at(i)};

        const std::size_t possible_weight {current_weight + current_item.weight};
        if (possible_weight <= capacity &&
            possible_weight * model_params.speed_to_weight_ratio < maximum_load_threshold) {
            current_profit += current_item.profit;
            current_weight += current_item.weight;
            knapsack_solution.at(current_item.node_index - 2) = true;
        }
    }

    solution.knapsack_value = current_profit;

    return knapsack_solution;
}

PackingPlan Model::solve_knapsack_dp(Solution& solution)
{
    const std::size_t items_count {items.size()};
    const std::size_t capacity {model_params.capacity};
    std::vector<std::vector<int>> dp_matrix(items_count + 1, std::vector<int>(capacity + 1));
    std::vector<bool> knapsack_solution(items_count);

    for (std::size_t i {0}; i <= items_count; i++) {
        for (std::size_t w {0}; w <= capacity; w++) {
            if (i == 0 || w == 0)
                dp_matrix.at(i).at(w) = 0;
            else if (items.at(i - 1).weight <= w)
                dp_matrix.at(i).at(w) = std::max(solution.penalized_items.at(i - 1).profit +
                    dp_matrix.at(i - 1).at(w - items.at(i - 1).weight), dp_matrix.at(i - 1).at(w));
            else
                dp_matrix.at(i).at(w) = dp_matrix.at(i - 1).at(w);
        }
    }

    int current_profit {dp_matrix.at(items_count).at(capacity)};

    solution.knapsack_value = current_profit;
    std::size_t current_weight {capacity};

    for (std::size_t i {items_count}; i > 0 && current_profit > 0; i--) {
        if (current_profit != dp_matrix.at(i - 1).at(current_weight)) {
            knapsack_solution.at(items.at(i - 1).index - 1) = true;

            current_profit = current_profit - solution.penalized_items.at(i - 1).profit;
            current_weight = current_weight - items.at(i - 1).weight;
        }
    }

    return knapsack_solution;
}

std::vector<Node> Model::process_invert_mutation(const std::vector<Node>& route)
{
    std::vector<Node> inverted_route {route};
    std::uniform_int_distribution<std::size_t> int_distribution(0, model_params.dimension);
    std::size_t first_inverse_point {int_distribution(rng)};
    std::size_t second_inverse_point {int_distribution(rng)};

    if (first_inverse_point > second_inverse_point)
        std::swap(first_inverse_point, second_inverse_point);

    std::reverse(
        inverted_route.begin() + first_inverse_point,
        inverted_route.begin() + second_inverse_point);

    return inverted_route;
}

std::vector<Item> Model::penalize_item_values(const std::vector<Node>& route)
{
    std::vector<Item> penalized_items {items};

    for (std::size_t i {0uz}; i < route.size(); i++) {
        const std::size_t current_node {route.at(i).index};

        if (current_node != 1) {
            Item& current_item {penalized_items.at(current_node - 2)};
            current_item.profit *= ln_factor(i + 2, route.size() + 1);
            current_item.ratio =
                static_cast<double>(current_item.profit) /
                static_cast<double>(current_item.weight);
        }
    }

    return penalized_items;
}
