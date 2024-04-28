#include "./include/vertex.h"
#include "./include/graph.h"
#include "./include/utility_operators.h"

#include <iostream>

int main()
{
    Graph graph(10);

    std::cout << "\nSize: " << graph.vertices.size();
    std::cout << "\nVertices: [" << graph.vertices << "]\n";

    // std::cout << "\nFind 3: " << graph.vertices.count(Vertex(3));
    // std::cout << "\nFind 10: " << graph.vertices.count(Vertex(10));
    // std::cout << "\nFind 100: " << graph.vertices.count(Vertex(100));

    graph.add_edge(1, 3);
    graph.add_edge(1, 2);
    graph.add_edge(1, 5);
    graph.add_edge(1, 3);

    for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
        std::cout << graph.vertices.at(i).print_neighbours();

    std::cout << "\n\n";

    graph.vertices.at(1).update_color(1);
    graph.vertices.at(3).update_color(2);

    std::cout << "\nVertices: [" << graph.vertices << "]\n";

    for (std::size_t i {0uz}; i < graph.vertices.size(); i++)
        std::cout << graph.vertices.at(i).print_neighbours();

    std::cout << "\n\n";

    return 0;
}
