FLAGS=-W -Wall -Werror -lGL -lGLU -lglut -lm -Ofast -ggdb
all: parallel nonparallel

parallel: prt pcrt

nonparallel: rt rc

prt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp
pcrt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp -DRAYTRACE_CACHE

rt: c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -DRAYTRACE_NONPARALLEL
rc: c/ray-cast.c c/ray-cast.h
	gcc $< $(FLAGS) -o bin/$@
clean:
	rm *