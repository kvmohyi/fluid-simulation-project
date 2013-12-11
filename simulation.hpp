#ifndef	_SIMULATION_HPP
#define _SIMULATION_HPP

#include <vector>

#include <glm/glm.hpp>

#include "geometry.hpp"

#define DEBUG true

using namespace std;
using namespace glm;

class FluidSimulation {
public:
	vector<vector<Particle> > gridCells;
	RigidBody cube;

	vec3 gravity; // by default, set this to (0, -9.8, 0)
	float timeStepSize; // the size of each timestep. the smaller, the more accurate.
	int numGrids; // the number of grids along each axis
	int numParticles; // the total number of particles in the simulation
	float worldSize; // the size of the world
	float gridSize; // the size of each side of each grid cell
	float localRadius; // the region of support for the kernel smoothing functions
	float viscosityConstant; // viscosity coefficient to be used in viscosity equation
	float gasConstant; // gas coefficient to be used in pressure equation
	float restDensity; // the rest density of the fluid. 1000 for water
	float particleMass; // the mass of each particle: volume * density / num_particles
	int testVersion; // 1 is a cube in the center, 2 is a cube dropping into water
	int dimensions; // set to 2 for 2D, 3 for 3D
	float tensionConstant;
	float tensionThreshold;

	int numIterations;

	void instantiateFromFile(string file);

	FluidSimulation();
	void elapseTimeGrid();
	vector<vector<Particle> >& particleList();
	int mapToIndex(Particle particle);
	int mapToIndex(Particle particle, int x_offset, int y_offset, int z_offset);
	int mapToIndex(int x, int y, int z);
    void drawWaterShape(int numParticles, float xStart, float yStart, float zStart, float xEnd, float yEnd, float zEnd);
    void drawTest(int dimension, int version);
    void printParams();
    double sphereRadius(Particle& particle);
};


#endif
