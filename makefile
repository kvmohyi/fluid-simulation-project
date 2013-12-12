all:
	g++ -fopenmp -I./glm main.cpp -I./FreeImage geometry.cpp equations.cpp simulation.cpp -o fluidsim -lglut -lGL -lGLU -lfreeimage