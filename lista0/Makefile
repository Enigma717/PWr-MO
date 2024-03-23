CPPFLAGS = -Wall -Wextra -pedantic -O3 -std=c++23

.PHONY = all clean

SRC  := $(wildcard src/*.cpp)

all:
	$(CXX) $(CPPFLAGS) -o sim main.cpp $(SRC)

clean:
	rm -f sim
