#include "../include/model.h"
#include "../include/graph.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <sstream>

Model::Model() : loader{*this}, graph(nullptr), solver{*this}, rng {rd()} {}

void Model::load_file(const std::string& file_path)
{
    loader.parse_instance(file_path);
}

void Model::create_graph()
{
    graph = std::make_unique<Graph>(model_params.vertices);
}

void Model::add_edge_to_graph(
    const std::size_t source_id,
    const std::size_t destination_id)
{
    graph->add_edge(source_id, destination_id);
}

std::size_t Model::calculate_max_degree() const
{
    std::size_t max_degree {0};

    for (const auto& vertex : graph->vertices) {
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
        << "\n|-> Maximum degree: " << model_params.max_degree;

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

std::size_t Model::evaluate_fitness(const Graph& solution)
{
    for (const auto* const colour : solution.colours)
        final_colours.insert(*colour);

    return check_colouring_corretness(solution) ? final_colours.size() : 999'999'999;
}

Graph Model::solve_random(Graph graph)
{
    graph.reset_colouring();
    solver.random_solution(graph);

    return graph;
}

Graph Model::solve_greedy(Graph graph)
{
    graph.reset_colouring();
    solver.greedy_solution(graph);

    return graph;
}
