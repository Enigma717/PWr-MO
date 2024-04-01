#include "./include/model.h"
#include "./include/enums/population_type.h"
#include "./include/utility_operators.h"

#include <chrono>
#include <iostream>

int main()
{
    Model model;

    model.load_file("./berlin52/berlin52_n51_uncorr_05.ttp");
    // model.load_file("./berlin52/berlin52_n51_uncorr-similar-weights_01.ttp");
    // model.load_file("./berlin52/berlin52_n51_bounded-strongly-corr_01.ttp");
    // model.load_file("./bier127/bier127_n126_uncorr_01.ttp");

    // model.load_file("./a280/a280_n279_uncorr_01.ttp");

    std::cout << "\n===========================\n\n"
        << model.print_model_parms() << "\n";

    model.create_weight_matrix();

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto best_member {model.genetic_solver.solve()};
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    std::cout << "\n========[ RESULTS ]========\n\n";
    std::cout << "|-> Best route: " << best_member.route << "\n";
    std::cout << "|-> Best solution fitness: " << best_member.fitness << "\n";
    std::cout << "|-> Knapsack value: " << best_member.knapsack_value << "\n";
    std::cout << "|-> Traveling time: " << best_member.knapsack_value - best_member.fitness << "\n";
    std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    std::cout << "\n===========================\n\n";

    return 0;
}
