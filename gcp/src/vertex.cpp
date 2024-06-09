#include "../include/vertex.h"
#include "../include/utility_operators.h"

#include <sstream>

Vertex::Vertex(const std::size_t id) : id{id} {};

std::string Vertex::print_neighbours() const
{
    std::stringstream stream;

    stream << "\nNeighbourhood(" << *this << "): [" << neighbours << "]";

    return stream.str();
}

std::string Vertex::print_indirect_neighbours() const
{
    std::stringstream stream;

    stream << "\nNeighbourhood(" << *this << "): [" << indirect_neighbours << "]";

    return stream.str();
}

std::size_t Vertex::get_id() const
{
    return id;
}

std::size_t Vertex::get_colour() const
{
    return colour;
}

std::size_t& Vertex::get_colour_ref()
{
    return colour;
}

const std::set<Vertex*>& Vertex::get_neighbours() const
{
    return neighbours;
}

const std::set<Vertex*>& Vertex::get_indirect_neighbours() const
{
    return indirect_neighbours;
}

void Vertex::update_neighbourship_with(Vertex& new_neighbour)
{
    neighbours.insert(&new_neighbour);
}

void Vertex::update_indirect_neighbourship_with(Vertex& new_neighbour)
{
    indirect_neighbours.insert(&new_neighbour);
}

void Vertex::update_colour(std::size_t new_colour)
{
    colour = new_colour;
}
