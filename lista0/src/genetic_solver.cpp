#include "../include/genetic_solver.h"

GeneticSolver::GeneticSolver(Model& model_ref) : model_ref{model_ref} {}

void GeneticSolver::initialize_population(
    PopulationType population_type, std::uint32_t population_size)
{
    population.resize(population_size);

    switch (population_type) {
        case PopulationType::random: random_initialization(); break;
        case PopulationType::neighbour: neighbour_initialization(); break;
    }
}

void GeneticSolver::random_initialization()
{
    for (size_t population_member {0}; population_member < population_size; population_member++) {
    }
}

void GeneticSolver::neighbour_initialization()
{
    for (size_t population_member {0}; population_member < population_size; population_member++) {
    }
}
