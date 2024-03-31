#pragma once

#include "./enums/population_type.h"
#include "./structs/node.h"
#include "./structs/member.h"

#include <string>
#include <ostream>

inline std::ostream& operator<<(std::ostream& stream, const ProblemType& problem_type)
{
    switch (problem_type) {
    case ProblemType::uncorrelated:
        stream << "UNCORRELATED";
        break;
    case ProblemType::bounded:
        stream << "BOUNDED";
        break;
    case ProblemType::similar:
        stream << "SIMILAR";
        break;
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Node>& solution)
{
    stream << "[";

    for (std::size_t i {0u}; i < solution.size(); i++) {
        stream << solution.at(i).index;

        if (i != solution.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const Member& member)
{
    stream << member.solution << ", fitness: " << member.fitness;

    return stream;
}

inline auto operator<=>(const Member& member1, const Member& member2)
{
    return member1.fitness <=> member2.fitness;
}
