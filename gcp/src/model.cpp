#include "model.hpp"
#include "graph.hpp"
#include "utility_operators.hpp"

#include <iostream>
#include <sstream>

Model::Model()
: loader {*this},
  base_graph {nullptr},
  solver {*this},
  genetic_solver {*this},
  ltgomea_solver {*this},
  p3_solver {*this},
  rng {rd()}
{
}

void Model::load_file(const std::string& file_path)
{
    loader.parse_instance(file_path);
}

void Model::create_base_graph()
{
    base_graph = std::make_unique<Graph>(model_params.vertices);
}

void Model::add_edge_to_base_graph(
    const std::size_t source_id,
    const std::size_t destination_id)
{
    base_graph->add_edge(source_id, destination_id);
}

std::size_t Model::calculate_max_degree() const
{
    std::size_t max_degree {0};

    for (const auto& vertex : base_graph->vertices) {
        const std::size_t degree {vertex.get_neighbours().size()};

        if (max_degree < degree) {
            max_degree = degree;
        }
    }

    return max_degree;
}

std::string Model::print_model_parms() const
{
    std::stringstream stream;
    stream << "-> Model parameters:"
        << "\n|-> Name: " << model_params.instance_name
        << "\n|-> Vertices: " << model_params.vertices
        << "\n|-> Edges: " << model_params.edges
        << "\n|-> Maximum degree: " << model_params.max_degree
        << "\n\\-> Optimum colouring: " << model_params.optimum
        << "\n";

    return stream.str();
}

bool Model::check_colouring_corretness(const Graph& solution) const
{
    for (const auto& vertex : solution.vertices) {
        for (const auto* neighbour : vertex.get_neighbours()) {
            if (neighbour->get_colour() == vertex.get_colour())
                return false;
        }
    }

    return true;
}

bool Model::check_no_colouring_collisions(
    const std::vector<std::size_t>& forbidden_colours, std::size_t vertex_colour) const
{
    const auto it {std::find_if(
        forbidden_colours.begin(),
        forbidden_colours.end(),
        [=](const std::size_t forbidden_colour) {
            return forbidden_colour == vertex_colour;
        }
    )};

    return it == forbidden_colours.end();
}

std::set<std::size_t> Model::get_used_colours(const Graph& solution) const
{
    std::set<std::size_t> used_colours;

    for (const auto* const colour : solution.colours)
        used_colours.insert(*colour);

    return used_colours;
}

std::vector<std::size_t> Model::get_forbidden_colours(const Vertex& vertex) const
{
    const auto& neighbours {vertex.get_neighbours()};
    const auto& indirect_neighbours {vertex.get_indirect_neighbours()};

    std::set<std::size_t> forbidden_colours_set;
    for (const auto* const neighbour : neighbours)
        forbidden_colours_set.insert(neighbour->get_colour());

    for (const auto* const neighbour : indirect_neighbours)
        forbidden_colours_set.insert(neighbour->get_colour());

    return {forbidden_colours_set.begin(), forbidden_colours_set.end()};
}

std::size_t Model::find_available_colour(const std::vector<std::size_t>& forbidden_colours) const
{
    std::size_t available_colour {1};

    if (forbidden_colours.size() == 0)
        return available_colour;

    if (forbidden_colours.front() > available_colour)
        return available_colour;

    for (std::size_t i {0}; i < forbidden_colours.size(); i++) {
        if (forbidden_colours.size() == 1)
            return forbidden_colours.at(0) == 1 ? 2 : 1;

        if (i == forbidden_colours.size() - 1)
            return forbidden_colours.at(i) + 1;

        if (forbidden_colours.at(i + 1) - forbidden_colours.at(i) > 1) {
            available_colour = (forbidden_colours.at(i) + 1);
            break;
        }
    }

    return available_colour;
}

Graph& Model::fix_colouring(Graph& solution) const
{
    for (auto& vertex : solution.vertices) {
        const auto forbidden_colours {get_forbidden_colours(vertex)};
        const std::size_t vertex_colour {vertex.get_colour()};

        if (check_no_colouring_collisions(forbidden_colours, vertex_colour))
            continue;

        const std::size_t available_colour {find_available_colour(forbidden_colours)};
        vertex.update_colour(available_colour);
    }

    return solution;
}

std::size_t Model::evaluate_fitness(Graph& solution) const
{
    return check_colouring_corretness(solution)
        ? get_used_colours(solution).size()
        : get_used_colours(fix_colouring(solution)).size();
}

void Model::mutate_random_vertex(Graph& graph)
{
    std::uniform_int_distribution<std::size_t> vertex_int_distribution(0, model_params.vertices - 1);
    std::uniform_int_distribution<std::size_t> colour_int_distribution(1, model_params.max_degree);

    const std::size_t random_vertex {vertex_int_distribution(rng)};
    const std::size_t random_colour {colour_int_distribution(rng)};

    graph.vertices.at(random_vertex).update_colour(random_colour);
}

bool Model::are_two_solutions_same(const Solution& sol1, const Solution& sol2)
{
    const auto& sol1_colours {sol1.graph.colours};
    const auto& sol2_colours {sol2.graph.colours};

    return std::equal(
        sol1_colours.begin(), sol1_colours.end(),
        sol2_colours.begin(), sol2_colours.end(),
        [](const std::size_t* col1, const std::size_t* col2){
            return *col1 == *col2;
        });
}

bool Model::check_for_equality_in_cluster(
    const Solution& sol1,
    const Solution& sol2,
    const std::vector<std::size_t> cluster)
{
    const auto& sol1_colours {sol1.graph.colours};
    const auto& sol2_colours {sol2.graph.colours};

    for (const auto variable : cluster)
        if (*sol1_colours.at(variable) == *sol2_colours.at(variable))
            return true;

    return false;
}

Graph Model::solve_random()
{
    Graph solution_graph {*base_graph};
    solver.random_solution(solution_graph);

    return solution_graph;
}

Graph Model::solve_random(Solution& solution)
{
    solution.graph.reset_colouring();
    solver.random_solution(solution.graph);

    return solution.graph;
}

Graph Model::solve_greedy()
{
    Graph solution_graph {*base_graph};
    solver.greedy_solution(solution_graph);

    return solution_graph;
}

Graph Model::solve_greedy(Solution& solution)
{
    solution.graph.reset_colouring();
    solver.greedy_solution(solution.graph);

    return solution.graph;
}

Graph Model::solve_simulated_annealing()
{
    Graph solution {*base_graph};
    solver.simulated_annealing_solution(solution);

    return solution;
}

void Model::solve_genetic()
{
    genetic_solver.solve();
}

void Model::solve_ltgomea()
{
    ltgomea_solver.solve();
}

void Model::solve_p3()
{
    p3_solver.solve();
}
