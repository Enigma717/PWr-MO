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

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Node>& route)
{
    stream << "[";

    for (std::size_t i {0u}; i < route.size(); i++) {
        stream << route.at(i).index;

        if (i != route.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Item>& items)
{
    stream << "[";

    for (std::size_t i {0uz}; i < items.size(); i++) {
        const auto& item {items.at(i)};
        stream << "\nItem(" << item.index << "): [node_idx: "
            << item.node_index << ", profit: "
            << item.profit << ", weight: "
            << item.weight << ", ratio: "
            << item.ratio << "]";

        if (i != items.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<bool>& packing_plan)
{
    stream << "[";

    for (std::size_t i {0u}; i < packing_plan.size(); i++) {
        stream << packing_plan.at(i);

        if (i != packing_plan.size() - 1)
            stream << ", ";
    }

    stream << "]";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const Member& member)
{
    stream << "Route: " << member.route
        << ", packing plan: " << member.packing_plan
        << ", fitness: " << member.fitness;

    return stream;
}

inline auto operator<=>(const Member& member1, const Member& member2)
{
    return member1.fitness <=> member2.fitness;
}

inline double operator+(const Member& member1, const Member& member2)
{
    return member1.fitness + member2.fitness;
}

inline double operator+(double sum, const Member& member)
{
    return sum + member.fitness;
}
