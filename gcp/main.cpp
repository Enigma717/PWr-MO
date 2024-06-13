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
    model.load_file("./instances/queen8_8.col");
    // model.load_file("./instances/inithx.i.1.col");
    // model.load_file("./instances/latin_square_10.col");

    std::cout << model.print_model_parms();
    std::cout << "\nSize: " << model.base_graph->vertices.size() << "\n";
    // std::cout << "\nVertices: [" << model.graph->vertices << "]\n";
    // std::cout << "\nColours: [" << model.graph->colours << "]\n";

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Graph solution {model.solve_random()};
    // Graph solution {model.solve_greedy()};
    Solution& solution {model.solve_genetic()};
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    // std::cout << "\nVertices: [" << model.graph->vertices << "]";
    // std::cout << "\nColours: [" << model.graph->colours << "]\n";

    // for (std::size_t i {0uz}; i < model.graph->vertices.size(); i++)
    //     std::cout << model.graph->vertices.at(i).print_neighbours();

    std::cout << "\n\n========[ RESULTS ]========\n\n";
    std::cout << "|-> Solution: " << solution << "\n";
    std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
    std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
    std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    std::cout << "\n===========================\n\n";

    return 0;
}
