#include "./include/model.h"
#include "./include/graph.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main()
{
    Model model;

    // model.load_file("./instances/bipartite.col");
    model.load_file("./instances/queen5_5.col");
    // model.load_file("./instances/inithx.i.1.col");
    // model.load_file("./instances/latin_square_10.col");

    std::cout << model.print_model_parms();
    std::cout << "\nSize: " << model.graph->vertices.size() << "\n";
    // std::cout << "\nVertices: [" << model.graph->vertices << "]\n";
    // std::cout << "\nColours: [" << model.graph->colours << "]\n";

    // for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
        // graph.vertices.at(i).update_colour(i);

    // for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
    //     std::cout << graph.vertices.at(i).print_neighbours();

    // std::cout << "\n\n";

    // std::cout << "-> BFS: ";
    // model.graph->BFS(0);

    // std::cout << "-> DFS: ";
    // model.graph->DFS(0);

    // std::cout << "-> Bipartite: ";
    // model.graph->is_bipartite();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // Graph solution {model.solve_random(*model.graph)};
    Graph solution {model.solve_greedy(*model.graph)};
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    // std::cout << "\nVertices: [" << model.graph->vertices << "]";
    // std::cout << "\nColours: [" << model.graph->colours << "]\n";

    // for (std::size_t i {0uz}; i < model.graph->vertices.size(); i++)
    //     std::cout << model.graph->vertices.at(i).print_neighbours();

    std::cout << "\n\n========[ RESULTS ]========\n\n";
    std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution) << "\n";
    std::cout << "|-> Final colours: [" << model.final_colours << "]\n";
    std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    std::cout << "\n===========================\n\n";

    return 0;
}
