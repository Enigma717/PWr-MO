#include "./include/model.hpp"
#include "./include/graph.hpp"
#include "./include/utility_operators.hpp"

#include <chrono>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "\nMissing tournament size!\n";

        return 1;
    }

    std::size_t subpopulations_limit {static_cast<std::size_t>(std::stoi(argv[1]))};

    Model model;

    std::vector<std::string> files {
        "anna",
        // "david",
        // "huck",
        // "queen5_5",
        // "queen6_6",
        // "queen7_7",
        // "queen8_8",
        // "myciel3",
        // "myciel4",
        // "myciel5",
        // "myciel6",
        // "myciel7",
        // "zeroin.i.1"
    };

    for (const auto& filename : files) {
        for (int l = 0; l < 1; l++) {
            std::stringstream path;
            path << "./instances/" << filename << ".col";

            model.load_file(path.str());
            std::cout << "\n\n========[ Instance ]========\n\n";
            std::cout << model.print_model_parms();

            model.ltgomea_solver.subpopulations_limit = subpopulations_limit;

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            model.solve_ltgomea();
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

            std::cout << "\n\n========[ Results ]========\n\n";
            std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
            std::cout << "\n===========================\n\n";
        }
    }

    // for (const auto& filename : files) {
    //     for (int l = 0; l < 1; l++) {
    //         std::stringstream path;
    //         path << "./instances/" << filename << ".col";

    //         model.load_file(path.str());

    //         double avg {0.0};

    //         std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //         Solution solution {model.solve_genetic(avg)};
    //         std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //         const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    //         std::cout << "\n\n========[ Results ]========\n\n";
    //         std::cout << "|-> Solution: " << solution << "\n";
    //         std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
    //         std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
    //         std::cout << "|-> Trounament size: [" << model.genetic_solver.tournament_size << "]\n";
    //         std::cout << "|-> Crossing probability: [" << model.genetic_solver.crossing_probability << "]\n";
    //         std::cout << "|-> Mutation probability: [" << model.genetic_solver.mutation_probability << "]\n";
    //         std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    //         std::cout << "\n===========================\n\n";
    //     }
    // }

    // for (const auto& filename : files) {
    //     std::stringstream path;
    //     path << "./instances/" << filename << ".col";

    //     model.load_file(path.str());

    //     std::cout << "\n\n========[ Instance ]========\n\n";
    //     std::cout << model.print_model_parms();

    //     std::ofstream results_file;
    //     std::stringstream results_path;
    //     results_path << "./csv/tuning/" << filename << "_subgroup" << tournament_size << ".csv";
    //     results_file.open(results_path.str(), std::ios_base::app);

    //     results_file << "cx; mx; best; avg; time\n";

    //     model.genetic_solver.tournament_size = tournament_size;

    //     for (int j = 0; j < 5; j++) {
    //         for (int k = 0; k < 5; k++) {
    //             const auto cxp {0.5 + (0.1 * j)};
    //             const auto mxp {0.05 + (0.05 * k)};

    //             model.genetic_solver.crossing_probability = cxp;
    //             model.genetic_solver.mutation_probability = mxp;

    //             for (int l = 0; l < 1; l++) {
    //                 model.genetic_solver.generation_number = 0;
    //                 model.genetic_solver.fitness_evaluations = 0;

    //                 double avg {0.0};

    //                 std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    //                 Solution solution {model.solve_genetic(avg)};
    //                 std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //                 const auto elapsed_time {(std::chrono::duration<double>(end - begin))};

    //                 std::cout << "\n\n========[ Results ]========\n\n";
    //                 std::cout << "|-> Solution: " << solution << "\n";
    //                 std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
    //                 std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
    //                 std::cout << "|-> Trounament size: [" << model.genetic_solver.tournament_size << "]\n";
    //                 std::cout << "|-> Crossing probability: [" << model.genetic_solver.crossing_probability << "]\n";
    //                 std::cout << "|-> Mutation probability: [" << model.genetic_solver.mutation_probability << "]\n";
    //                 std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    //                 std::cout << "\n===========================\n\n";

    //                 results_file << cxp << "; "
    //                     << mxp << "; "
    //                     << model.evaluate_fitness(solution.graph) << "; "
    //                     << avg << "; "
    //                     << elapsed_time << "\n";
    //             }
    //         }
    //     }

    //     results_file.close();
    // }

    // std::cout << "\n\n========[ Results ]========\n\n";
    // std::cout << "|-> Solution: " << solution << "\n";
    // std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution.graph) << "\n";
    // std::cout << "|-> Final colours: [" << model.get_used_colours(solution.graph) << "]\n";
    // std::cout << "|-> Time elapsed: " << elapsed_time << "\n";
    // std::cout << "\n===========================\n\n";

    // for (const auto& filename : files) {
    //     for (int l = 0; l < 1; l++) {
    //         std::stringstream path;
    //         path << "./instances/" << filename << ".col";

    //         model.load_file(path.str());

    //         model.genetic_solver.generation_number = 0;
    //         model.genetic_solver.fitness_evaluations = 0;

    //         Graph solution {model.solve_random()};

    //         std::cout << "|-> Solution fitness: " << model.evaluate_fitness(solution) << "\n";
    //     }
    // }

    // for (const auto& filename : files) {
    //     for (int l = 0; l < 1; l++) {
    //         std::stringstream path;
    //         path << "./instances/" << filename << ".col";

    //         model.load_file(path.str());
    //         std::cout << "\n\n========[ Instance ]========\n\n";
    //         std::cout << model.print_model_parms();

    //         Graph solution {model.solve_greedy()};

    //         std::cout << "|-> Solution fitness GREEEEDY: " << model.evaluate_fitness(solution) << "\n";
    //     }
    // }

    return 0;
}
