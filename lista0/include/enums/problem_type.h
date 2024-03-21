#pragma once

#include <string>
#include <ostream>

enum class ProblemType {
    uncorrelated,
    bounded,
    similar
};

inline std::ostream& operator<<(std::ostream& stream, ProblemType& problem_type)
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
