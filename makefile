all:
	g++ -I./glm main.cpp simulation.cpp equations.cpp geometry.cpp -o fluidsim
testgen:
	g++ testgenerator.cpp -o testgen
