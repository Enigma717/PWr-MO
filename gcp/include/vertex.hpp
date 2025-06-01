#pragma once

#include <cstdint>
#include <set>
#include <string>

class Vertex {
public:
    Vertex() = delete;
    Vertex(const std::size_t);

    bool operator==(const Vertex&) const = default;

    std::string print_neighbours() const;
    std::string print_indirect_neighbours() const;
    std::size_t get_id() const;
    std::size_t get_colour() const;
    std::size_t& get_colour_ref();
    const std::set<Vertex*>& get_neighbours() const;
    const std::set<Vertex*>& get_indirect_neighbours() const;

    void update_neighbourship_with(Vertex& new_neighbour);
    void update_indirect_neighbourship_with(Vertex& new_neighbour);
    void update_colour(std::size_t new_colour);

// private:
    std::size_t id {0};
    std::size_t colour {0};
    std::set<Vertex*> neighbours;
    std::set<Vertex*> indirect_neighbours;
};
