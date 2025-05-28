#include "subpopulation.hpp"
#include "model.hpp"
#include "utility_operators.hpp"

#include <iostream>

Subpopulation::Subpopulation(std::size_t subpopulation_size, Model& model_ref)
: subpopulation_size {subpopulation_size},
  model_ref {model_ref},
  best_solution {nullptr},
  worst_solution {nullptr}
{
    individuals.resize(subpopulation_size);

    for (auto& solution : individuals)
        solution = create_new_solution(model_ref.solve_random());
}

void Subpopulation::print_individuals() const
{
    for (std::size_t i {0uz}; i < individuals.size(); i++)
        std::cout << "Subpopulation size: " << subpopulation_size << " | Solution " << i
            << " (" << &individuals.at(i) << "): " << individuals.at(i) << "\n";
}

void Subpopulation::print_info() const
{
    std::cout << "Subpopulation size: " << subpopulation_size
        << " | Best solution: [" << best_solution << "], fitness: " << best_solution->fitness
        << " | Worst solution: [" << worst_solution << "], fitness: " << worst_solution->fitness
        << " | Average fitness: " << avg_fitness
        << " | Iterations done: " << iterations_done << "\n";
}

void Subpopulation::run_iteration()
{
    for (auto& solution : individuals) {
        model_ref.mutate_random_vertex(solution.graph);
        solution.fitness = fitness_evaluation(solution);
    }

    auto new_best_solution {&*std::min_element(individuals.begin(), individuals.end())};
    auto new_worst_solution {&*std::max_element(individuals.begin(), individuals.end())};
    double new_avg_fitness {
        static_cast<double>(std::accumulate(individuals.begin(), individuals.end(), 0.0)) /
        static_cast<double>(subpopulation_size)};

    if (iterations_done == 0) {
        best_solution = new_best_solution;
        worst_solution = new_worst_solution;
    }

    if (new_best_solution->fitness <= best_solution->fitness)
        best_solution = new_best_solution;

    if (new_worst_solution->fitness >= worst_solution->fitness)
        worst_solution = new_worst_solution;

    avg_fitness = new_avg_fitness;

    iterations_done++;
}

Solution Subpopulation::create_new_solution(Graph&& graph)
{
    Solution solution(std::move(graph));
    solution.fitness = fitness_evaluation(solution);

    return solution;
}

double Subpopulation::fitness_evaluation(Solution& solution)
{
    fitness_evaluations++;
    return model_ref.evaluate_fitness(solution.graph);
}
