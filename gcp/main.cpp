#include "./include/model.h"
#include "./include/graph.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main()
{
    Model model;

    // model.load_file("./instances/bipartite.col");
    // model.load_file("./instances/queen5_5.col");
    // model.load_file("./instances/queen6_6.col");
    model.load_file("./instances/queen7_7.col");
    // model.load_file("./instances/queen8_8.col");
    // model.load_file("./instances/inithx.i.1.col");
    // model.load_file("./instances/latin_square_10.col");

    std::cout << "\n\n========[ Instance ]========\n\n";
    std::cout << model.print_model_parms();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Graph solution {model.solve_random()};
    // Graph solution {model.solve_greedy()};
    Graph solution {model.solve_simulated_annealing()};
    // Graph& solution {model.solve_genetic().graph};
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    std::cout << "\n\n========[ Results ]========\n\n";
    std::cout << "|-> Solution: " << solution << "\n";
    std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution) << "\n";
    std::cout << "|-> Final colours: [" << model.get_used_colours(solution) << "]\n";
    std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    std::cout << "\n===========================\n\n";

    return 0;
}
