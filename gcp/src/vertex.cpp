#include "../include/vertex.h"
#include "../include/utility_operators.h"

#include <sstream>

Vertex::Vertex() : id{next_id++} {};

std::string Vertex::print_neighbours() const
{
    std::stringstream stream;

    stream << "\nNeighbours(" << *this << "): [" << neighbours << "]";

    return stream.str();
}

std::size_t Vertex::get_id() const
{
    return id;
}

std::size_t Vertex::get_color() const
{
    return color;
}

void Vertex::update_neighbourship(Vertex& new_neighbour)
{
    neighbours.insert(&new_neighbour);
}

void Vertex::update_color(std::size_t new_color)
{
    color = new_color;
}
