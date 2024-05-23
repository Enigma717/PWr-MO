#include "./include/model.h"
#include "./include/graph.h"
#include "./include/utility_operators.h"

#include <iostream>

int main()
{
    Model model;

    // model.load_file("./instances/inithx.i.1.col");
    // model.load_file("./instances/queen5_5.col");
    model.load_file("./instances/bipartite.col");

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

    model.solve_random();

    std::cout << "\nVertices: [" << model.graph->vertices << "]";
    std::cout << "\nColours: [" << model.graph->colours << "]\n";

    for (std::size_t i {0uz}; i < model.graph->vertices.size(); i++)
        std::cout << model.graph->vertices.at(i).print_neighbours();

    std::cout << "\n\nFITNESS EVALUATION: " << model.evaluate_fitness() << "\n\n";
    return 0;
}
