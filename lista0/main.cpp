#include "./include/model.h"
#include "./include/enums/population_type.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main()
{
    Model model;

    model.load_file("./berlin52/berlin52_n51_uncorr_01.ttp");

    std::cout << "\n==========================\n\n"
        << model.print_model_parms() << "\n";

    // std::cout << "\n==========================\n\n"
    //     << "|-> Starting solution: " << model.print_nodes() << "\n";

    model.create_weight_matrix();
    // std::cout << "\n==========================\n\n"
    //     << "|-> Weights: " << model.print_weight_matrix() << "\n";


    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // // model.k_random_solution(10000);
    // // const auto solution {model.nearest_neighbour(20)};
    // const auto solution {model.extended_nearest_neighbour()};
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    // const auto distance {model.objective_function(solution)};
    // const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    model.genetic_solver.solve();

    // std::vector<Node> s1 = {
    //     {2, 0.0, 0.0}, {9, 0.0, 0.0}, {4, 0.0, 0.0},
    //     {8, 0.0, 0.0}, {1, 0.0, 0.0}, {6, 0.0, 0.0},
    //     {7, 0.0, 0.0}, {3, 0.0, 0.0}, {5, 0.0, 0.0},
    //     };
    // std::vector<Node> s2 = {
    //     {4, 0.0, 0.0}, {7, 0.0, 0.0}, {2, 0.0, 0.0},
    //     {3, 0.0, 0.0}, {9, 0.0, 0.0}, {5, 0.0, 0.0},
    //     {1, 0.0, 0.0}, {8, 0.0, 0.0}, {6, 0.0, 0.0},
    //     };

    // Member p1 {s1, 0.0};
    // Member p2 {s2, 0.0};
    // model.genetic_solver.order_crossover(p1, p2);
    // model.genetic_solver.order_crossover(p2, p1);

    // std::cout << "|-> Best solution: " << solution << "\n";
    // std::cout << "|-> Distance: " << distance << "\n";
    // std::cout << "|-> PRD: " << model.prd(distance) << "%\n";
    // std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    // std::cout << "\n==========================\n\n";

    return 0;
}
