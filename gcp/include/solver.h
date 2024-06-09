#pragma once

#include <memory>
#include <random>

class Model;
class Graph;
class Vertex;

class Solver {
public:
    Solver() = delete;
    Solver(Model& model_ref);

    void random_solution(Graph& graph) const;
    void greedy_solution(Graph& graph) const;

private:
    void colour_randomly(Vertex& vertex) const;
    std::size_t find_first_avail_colour(Vertex& vertex) const;

    Model& model_ref;
};
