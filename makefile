FLAGS=-W -Wall -Werror -lGL -lGLU -lglut -lm -O3
all: parallel nonparallel

parallel: prt

nonparallel: rt rc

prt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o $@ -fopenmp

rt: c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o $@ -DRAYTRACE_NONPARALLEL
rc: c/ray-cast.c c/ray-cast.h
	gcc $< $(FLAGS) -o $@
clean:
	rm *