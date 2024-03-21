#include "../include/genetic_solver.h"
#include "../include/model.h"

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
    for (auto& population_member : population) {
        population_member.solution = model_ref.k_random_solution(100);
        population_member.objective = model_ref.objective_function(population_member.solution);
    }
}

void GeneticSolver::neighbour_initialization()
{
    for (size_t member_position {0}; member_position < population.size(); member_position++) {
        auto& population_member {population.at(member_position)};
        const std::uint16_t starting_node_index {
            static_cast<std::uint16_t>((member_position % model_ref.model_params.dimension) + 1)};

        population_member.solution = model_ref.nearest_neighbour(starting_node_index);
        population_member.objective = model_ref.objective_function(population_member.solution);
    }
}
