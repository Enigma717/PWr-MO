#pragma once

#include "structs/model_params.h"
#include "loader.h"

#include <memory>
#include <random>

class Graph;

class Model {
public:
    Model();

    void load_file(const std::string& file_path);
    void create_graph();
    void add_edge_to_graph(
        const std::size_t source_id,
        const std::size_t destination_id);

    std::string print_model_parms() const;

private:
    ModelParams model_params;
    Loader loader;
    std::unique_ptr<Graph> graph;

    std::random_device rd;
    std::mt19937 rng;

    friend class Loader;
};
