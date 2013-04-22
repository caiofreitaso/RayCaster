FLAGS=-W -Wall -Werror -lGL -lGLU -lglut -lm
all: rt rc

rt: c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o $@
rc: c/ray-cast.c c/ray-cast.h
	gcc $< $(FLAGS) -o $@
clean:
	rm *