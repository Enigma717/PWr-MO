#pragma once

#include "./vertex.h"

#include <ostream>
#include <vector>

inline std::ostream& operator<<(std::ostream& stream, const Vertex& vertex)
{
    stream << "V(id: " << vertex.get_id() << ", color: " << vertex.get_color() << ")";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Vertex>& vertices)
{
    for (std::size_t i {0uz}; i < vertices.size(); i++) {
        const auto& vertex {vertices.at(i)};

        stream << vertex;

        if (i != vertices.size() - 1)
            stream << ", ";
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::set<Vertex*>& neighbours)
{
    for (auto it {neighbours.begin()}; it != neighbours.end(); it++) {

        stream << **it;

        if (it != std::prev(neighbours.end(), 1))
            stream << ", ";
    }

    return stream;
}

inline bool operator<(const Vertex& lhs, const Vertex& rhs)
{
    return lhs.get_id() < rhs.get_id();
}
