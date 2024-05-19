#include "./include/model.h"
#include "./include/utility_operators.h"

#include <iostream>

int main()
{
    Model model;

    // model.load_file("./instances/inithx.i.1.col");
    model.load_file("./instances/queen5_5.col");
    std::cout << model.print_model_parms() << "\n";

    std::cout << "\nSize: " << model.graph->vertices.size();
    std::cout << "\nVertices: [" << model.graph->vertices << "]";
    std::cout << "\nColours: [" << model.graph->colours << "]\n";

    // for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
        // graph.vertices.at(i).update_colour(i);

    // for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
    //     std::cout << graph.vertices.at(i).print_neighbours();

    // std::cout << "\n\n";

    // std::cout << "-> BFS: ";
    // graph.BFS(0);

    // std::cout << "-> DFS: ";
    // graph.DFS(0);

    // std::cout << "-> Bipartite: ";
    // graph.is_bipartite();

    return 0;
}
