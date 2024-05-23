#pragma once

#include "./vertex.h"

#include <ostream>
#include <vector>

inline std::ostream& operator<<(std::ostream& stream, const Vertex& vertex)
{
    stream << "(i:" << vertex.get_id() + 1 << ", c:" << vertex.get_colour() << ")";

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::set<std::size_t>& colours)
{
    const std::size_t size {colours.size()};
    std::size_t i {0};

    for (const auto colour : colours) {
        stream << colour;

        if (i != size - 1)
            stream << ", ";

        i++;
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<std::size_t*>& colours)
{
    const std::size_t size {colours.size()};

    for (std::size_t i {0uz}; i < size; i++) {
        const auto& colour {colours.at(i)};

        stream << *colour;

        if (i != size - 1)
            stream << ", ";
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Vertex>& vertices)
{
    const std::size_t size {vertices.size()};

    for (std::size_t i {0uz}; i < size; i++) {
        const auto& vertex {vertices.at(i)};

        stream << vertex;

        if (i != size - 1)
            stream << ", ";
    }

    return stream;
}

inline std::ostream& operator<<(std::ostream& stream, const std::vector<Vertex*>& vertices)
{
    const std::size_t size {vertices.size()};

    for (std::size_t i {0uz}; i < size; i++) {
        const auto& vertex {vertices.at(i)};

        stream << *vertex;

        if (i != size - 1)
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
