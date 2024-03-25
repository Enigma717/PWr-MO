#include "../include/loader.h"
#include "../include/model.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>

namespace
{
    constexpr std::size_t header_size {8};
    constexpr std::size_t coord_section_start {9};

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

    if(!input_stream.is_open()) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    while (getline(input_stream, read_line)) {
        if (currently_read_line < header_size)
            dispatch_header_parser(read_line);
        else if (currently_read_line > coord_section_start && currently_read_line < coord_section_end)
            push_node_into_vector(read_line);
        else if (currently_read_line > coord_section_end)
            push_item_into_vector(read_line);

        currently_read_line++;
    }

}

void Loader::dispatch_header_parser(const std::string& read_line)
{
    const auto delimiter_pos {read_line.find('\t')};
    const auto value {read_line.substr(delimiter_pos + 1, std::string::npos)};

    switch (currently_read_line) {
    case 0:
        model_ref.model_params.instance_name = value;
        break;
    case 1:
        model_ref.model_params.problem_type = decide_problem_type(read_line);
        break;
    case 2:
        model_ref.model_params.dimension = std::stoi(value);
        coord_section_end = coord_section_start + model_ref.model_params.dimension + 1;
        model_ref.nodes.reserve(model_ref.model_params.dimension);
        model_ref.weights.resize(model_ref.model_params.dimension);

        for (auto i {0}; i < model_ref.model_params.dimension; i++)
        {
            model_ref.weights.at(i).resize(model_ref.model_params.dimension);
        }

        break;
    case 3:
        model_ref.model_params.num_of_items = std::stoi(value);
        model_ref.items.reserve(model_ref.model_params.num_of_items);
        break;
    case 4:
        model_ref.model_params.capacity = std::stoi(value);
        break;
    case 5:
        model_ref.model_params.min_speed = std::stod(value);
        break;
    case 6:
        model_ref.model_params.max_speed = std::stod(value);
        break;
    case 7:
        model_ref.model_params.renting_ratio = std::stod(value);
        break;
    default: break;
    }
}

void Loader::push_node_into_vector(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, '\t')};

    Node node;
    node.index = std::stoi(tokens.at(0));
    node.x_coord = std::stod(tokens.at(1));
    node.y_coord = std::stod(tokens.at(2));

    model_ref.nodes.push_back(node);
}

void Loader::push_item_into_vector(const std::string& read_line)
{
    const auto tokens {tokenize(read_line, '\t')};

    Item item;
    item.index = std::stoi(tokens.at(0));
    item.profit = std::stoi(tokens.at(1));
    item.weight = std::stoi(tokens.at(2));
    item.node_index = std::stoi(tokens.at(3));

    model_ref.items.push_back(item);
}

ProblemType Loader::decide_problem_type(const std::string& read_line)
{
    const auto delimiter_pos {read_line.find(':')};
    const auto read_problem_type {read_line.substr(delimiter_pos + 2, std::string::npos)};

    if (read_problem_type == "uncorrelated\r")
        return ProblemType::uncorrelated;
    else if (read_problem_type == "uncorrelated, similar weights\r")
        return ProblemType::similar;
    else
        return ProblemType::bounded;
}
