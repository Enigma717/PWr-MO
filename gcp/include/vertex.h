#pragma once

#include <cstdint>
#include <set>
#include <string>

class Vertex {
public:
    Vertex();

    std::string print_neighbours() const;
    std::size_t get_id() const;
    std::size_t get_color() const;

    void update_neighbourship(Vertex& new_neighbour);
    void update_color(std::size_t new_color);

private:
    static inline std::size_t next_id {0};
    std::size_t id {0};
    std::size_t color {0};
    std::set<Vertex*> neighbours;
};
