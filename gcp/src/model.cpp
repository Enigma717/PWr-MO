#include "../include/model.h"

#include "../include/graph.h"
#include "../include/utility_operators.h"

#include <iostream>
#include <sstream>

Model::Model() : loader{*this}, graph(nullptr), rng {rd()} {}

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

std::string Model::print_model_parms() const
{
    std::stringstream stream;
    stream << "-> Model parameters:"
        << "\n|-> Name: " << model_params.instance_name
        << "\n|-> Vertices: " << model_params.vertices
        << "\n|-> Edges: " << model_params.edges;

    return stream.str();
}
