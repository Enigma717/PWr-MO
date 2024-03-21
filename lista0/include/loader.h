#pragma once

#include "./enums/problem_type.h"

class Model;

class Loader {
public:
    Loader() = delete;
    Loader(Model& model_ref);

    void parse_instance(const std::string& file_path);

private:
    Model& model_ref;

    int currently_read_line {0};
    int coord_section_end {9};


    void dispatch_header_parser(const std::string& read_line);
    void push_node_into_vector(const std::string& read_line);
    void push_item_into_vector(const std::string& read_line);
    ProblemType decide_problem_type(const std::string& read_line);
};
