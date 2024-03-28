#include "./include/model.h"
#include "./include/enums/population_type.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main()
{
    Model model;

    // model.load_file("./berlin52/berlin52_n51_uncorr_01.ttp");
    model.load_file("./a280/a280_n279_bounded-strongly-corr_01.ttp");

    std::cout << "\n==========================\n\n"
        << model.print_model_parms() << "\n";

    // std::cout << "\n==========================\n\n"
    //     << "|-> Starting solution: " << model.print_nodes() << "\n";

    model.create_weight_matrix();
    // std::cout << "\n==========================\n\n"
    //     << "|-> Weights: " << model.print_weight_matrix() << "\n";


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // model.k_random_solution(10000);
    // const auto solution {model.nearest_neighbour(2)};
    const auto solution {model.genetic_solver.solve()};
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    const auto distance {model.objective_function(solution)};
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    std::cout << "|-> Best solution: " << solution << "\n";
    std::cout << "|-> Distance: " << distance << "\n";
    std::cout << "|-> PRD: " << model.prd(distance) << "%\n";
    std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    std::cout << "\n==========================\n\n";

    return 0;
}
