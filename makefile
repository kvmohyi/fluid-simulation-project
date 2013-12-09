all:
	g++ -I./glm main.cpp simulation.cpp equations.cpp geometry.cpp -o fluidsim -lglut -lGL
