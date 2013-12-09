#ifndef	_SIMULATION_HPP
#define _SIMULATION_HPP

#include <vector>

#include <glm/glm.hpp>

#include "geometry.hpp"

using namespace std;
using namespace glm;

class FluidSimulation {
private:
	vector<vector<Particle> > gridCells;
	vector<RigidBody> rigidBodies;

	vec3 gravity; // by default, set this to (0, 0, -9.8)
	float timeStepSize; // the size of each timestep. the smaller, the more accurate.
	int numGrids; // the number of grids along each axis
	int numParticles; // the total number of particles in the simulation
	float worldSize; // the size of the world
	float gridSize; // the size of each side of each grid cell
	float localRadius; // the region of support for the kernel smoothing functions
	float volume; // the total volume of fluid
	float viscosityConstant; // viscosity coefficient to be used in viscosity equation
	float gasConstant; // gas coefficient to be used in pressure equation
	float restDensity; // the rest density of the fluid. 1000 for water
	float particleMass; // the mass of each particle: volume * density / num_particles

	void instantiateFromFile(string file);

public:
	FluidSimulation(string file);
	void elapseTimeGrid();
	vector<Particle>& particleList();
	int mapToIndex(Particle particle);
	int mapToIndex(Particle particle, int x_offset, int y_offset, int z_offset);
	int mapToIndex(int x, int y, int z);
};

#endif