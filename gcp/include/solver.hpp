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
    void simulated_annealing_solution(Graph& graph) const;

    double sa_initial_temperature {50'000.0};
    double sa_cooling_rate {0.995};


private:
    bool check_reached_optimum(const std::size_t fitness) const;
    bool check_reached_iteration_limit(const std::size_t iteration) const;
    void colour_vertex_randomly(Vertex& vertex) const;
    void colour_vertex_greedily(Vertex& vertex) const;

    Model& model_ref;
};
