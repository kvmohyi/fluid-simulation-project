all:
	g++ -I./glm main.cpp geometry.cpp equations.cpp simulation.cpp -o fluidsim -lglut -lGL -lGLU