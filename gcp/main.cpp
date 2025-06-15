#include "./include/model.hpp"
#include "./include/graph.hpp"
#include "./include/utility_operators.hpp"

#include <chrono>
#include <string>
#include <sstream>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "\nMissing needed parameters!\n";
        std::cerr << "\nCorrect usage: ./gcp <instance_name> <optimizer>\n";
        std::cerr << "\n<optimizer> options: 1 - GA-PX, 2 - LT-GOMEA, 3 - LT-GOMEA-PX, 4 - P3, 5 - P3-PX\n";

        return 1;
    }

    std::string instance_name {argv[1]};
    std::size_t optimizer_choice {static_cast<std::size_t>(std::stoi(argv[2]))};


    Model model;

    std::stringstream path;
    path << "./instances/" << instance_name << ".col";

    model.load_file(path.str());
    std::cout << "\n\n========[ Instance ]========\n\n";
    std::cout << model.print_model_parms();


    switch (optimizer_choice) {
    case 1:
        if (argc < 4) {
            std::cerr << "\nMissing GA parameters!\n";
            std::cerr << "\nCorrect usage: ./gcp <instance_name> 1 <crossing_prob> <mutation_prob>\n";

            return 1;
        }

        double crossing_prob {static_cast<double>(std::stod(argv[3]))};
        double mutation_prob {static_cast<double>(std::stod(argv[4]))};

        if (crossing_prob < 0.0 || crossing_prob > 1.0 || mutation_prob < 0.0 || mutation_prob > 1.0) {
            std::cerr << "\nIncorrect probabilities!\n";

            return 1;
        }

        model.genetic_solver.crossing_probability = crossing_prob;
        model.genetic_solver.mutation_probability = mutation_prob;

        model.solve_genetic();

        break;
    case 2:
        model.ltgomea_solver.crossover_type = CrossoverType::optimal_mixing;
        break;
    case 3:
        model.solve_ltgomea();
        break;
    case 4:
        model.p3_solver.crossover_type = CrossoverType::optimal_mixing;
        break;
    case 5:
        model.solve_p3();
        break;
    default:
        break;
    }

    return 0;
}
