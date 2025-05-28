#include "loader.hpp"
#include "model.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>

namespace
{
    std::vector<std::string> tokenize(const std::string& string, const char delimiter)
    {
        std::vector<std::string> tokens;
        std::stringstream temp(string);
        std::string token;

        while (getline(temp, token, delimiter)) {
            tokens.push_back(token);
        }

        return tokens;
    }
}

Loader::Loader(Model& model_ref) : model_ref{model_ref} {}

void Loader::parse_instance(const std::string& file_path)
{
    std::string read_line;
    std::ifstream input_stream;
    input_stream.open(file_path);

    std::cout << "\nFILE PATH: " << file_path;

    if(!input_stream.is_open()) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    reading_first_line = true;

    while (getline(input_stream, read_line)) {
        const char& line_type {read_line.at(0)};

        switch (line_type) {
        case 'c': parse_header(read_line); break;
        case 'p': parse_metadata(read_line); break;
        case 'o': parse_optimum(read_line); break;
        case 'e': parse_edge(read_line); break;
        }
    }

    model_ref.model_params.max_degree = model_ref.calculate_max_degree();
}

void Loader::parse_header(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, ' ')};

    if (reading_first_line) {
        model_ref.model_params.instance_name = tokens.back();
        reading_first_line = false;
    }
}

void Loader::parse_metadata(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, ' ')};

    model_ref.model_params.vertices = std::stoi(tokens.at(2));
    model_ref.model_params.edges = std::stoi(tokens.at(3));

    model_ref.create_base_graph();
}

void Loader::parse_optimum(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, ' ')};

    model_ref.model_params.optimum = std::stoi(tokens.at(1));
}

void Loader::parse_edge(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, ' ')};
    const auto source_vertex_id {std::stoi(tokens.at(1)) - 1};
    const auto destination_vertex_id {std::stoi(tokens.at(2)) - 1};

    model_ref.add_edge_to_base_graph(source_vertex_id, destination_vertex_id);
}
