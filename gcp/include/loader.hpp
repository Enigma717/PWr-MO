#pragma once

#include <cstdint>
#include <string>

class Model;

class Loader {
public:
    Loader() = delete;
    Loader(Model& model_ref);

    void parse_instance(const std::string& file_path);

private:
    void parse_header(const std::string& read_line);
    void parse_metadata(const std::string& read_line);
    void parse_optimum(const std::string& read_line);
    void parse_edge(const std::string& read_line);

    Model& model_ref;

    bool reading_first_line {true};
};
