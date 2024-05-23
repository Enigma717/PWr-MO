#pragma once

#include <memory>
#include <random>

class Model;
class Vertex;

class Solver {
public:
    Solver() = delete;
    Solver(Model& model_ref);

    void random_solution() const;

private:
    void draw_colour(Vertex& vertex) const;

    Model& model_ref;
};
