FLAGS=-W -Wall -Werror -lGL -lGLU -lglut -lm -Ofast -ggdb
all: parallel nonparallel

parallel: prt pcrt bprt

nonparallel: rt crt rc brt

prt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp
pcrt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp -DRAYTRACE_CACHE
bprt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp -DRAYTRACE_BRUTALMODE

rt: c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -DRAYTRACE_NONPARALLEL
crt: c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -DRAYTRACE_NONPARALLEL -DRAYTRACE_CACHE
brt:c++/raytrace.cpp c++/raytrace.hpp
	g++ $< $(FLAGS) -Wno-reorder -o bin/$@ -fopenmp -DRAYTRACE_BRUTALMODE -DRAYTRACE_NONPARALLEL

rc: c/ray-cast.c c/ray-cast.h
	gcc $< $(FLAGS) -o bin/$@
clean:
	rm bin/*