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
    void colour_vertex_randomly(Vertex& vertex) const;
    void colour_vertex_greedily(Vertex& vertex) const;

    Model& model_ref;
};
