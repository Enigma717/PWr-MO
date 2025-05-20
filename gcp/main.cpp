#include "./include/model.h"
#include "./include/graph.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main(int argc, char* argv[])
{
    Model model;

    // model.load_file("./instances/queen5_5.col");
    // model.load_file("./instances/queen6_6.col");
    // model.load_file("./instances/queen7_7.col");
    // model.load_file("./instances/queen8_8.col");
    // model.load_file("./instances/queen9_9.col");
    // model.load_file("./instances/queen10_10.col");

    // std::cout << "\n\n========[ Instance ]========\n\n";
    // std::cout << model.print_model_parms();

    // Graph solution {model.solve_random()};
    // Graph solution {model.solve_greedy()};
    // Graph solution {model.solve_simulated_annealing()};

    std::vector<std::string> files {
        // "queen5_5"
        // "queen6_6"
        // "queen7_7",
        "queen8_8"
        // "myciel5"
        // "myciel6",
        // "zeroin.i.1"
    };

    for (const auto& filename : files) {
        std::stringstream path;
        path << "./instances/" << filename << ".col";

        model.load_file(path.str());

        std::cout << "\n\n========[ Instance ]========\n\n";
        std::cout << model.print_model_parms();
    }

    for (const auto& filename : files) {
        for (int l = 0; l < 1; l++) {
            std::stringstream path;
            path << "./instances/" << filename << ".col";

            model.load_file(path.str());

            std::cout << "\n\n========[ Instance ]========\n\n";
            std::cout << model.print_model_parms();

            model.genetic_solver.generation_number = 0;
            model.genetic_solver.fitness_evaluations = 0;

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            Solution solution {model.solve_genetic()};
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

            std::cout << "\n\n========[ Results ]========\n\n";
            std::cout << "|-> Solution: " << solution << "\n";
            std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
            std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
            std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
            std::cout << "\n===========================\n\n";
        }
    }

    // std::cout << "\n\n========[ Results ]========\n\n";
    // std::cout << "|-> Solution: " << solution << "\n";
    // std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
    // std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
    // std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    // std::cout << "\n===========================\n\n";

    // for (const auto& filename : files) {
        // for (int l = 0; l < 10; l++) {
            // std::stringstream path;
            // path << "./instances/" << filename << ".col";
//
            // model.load_file(path.str());
//
            // model.genetic_solver.generation_number = 0;
            // model.genetic_solver.fitness_evaluations = 0;
//
            // Graph solution {model.solve_random()};
//
            // std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution) << "\n";
        // }
    // }
//
        // for (const auto& filename : files) {
        // for (int l = 0; l < 1; l++) {
            // std::stringstream path;
            // path << "./instances/" << filename << ".col";
//
            // model.load_file(path.str());
//
            // model.genetic_solver.generation_number = 0;
            // model.genetic_solver.fitness_evaluations = 0;
//
            // Graph solution {model.solve_greedy()};
//
            // std::cout << "|-> Solution fitness GREEEEDY: " << model.evaluate_fitness(solution) << "\n";
        // }
    // }

    return 0;
}
