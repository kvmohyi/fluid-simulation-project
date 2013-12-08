#ifndef	_SIMULATION_HPP
#define _SIMULATION_HPP

#include <vector>

#include <glm/glm.hpp>

#include "geometry.hpp"

using namespace std;
using namespace glm;

class FluidSimulation {
private:
	vector<Vector<Particle> > buckets;
	vector<RigidBody> rigidBodies;

	vec3 gravity;
	float timeStepSize;
	int numGrids;
	int numParticles;
	float gridSize;
	float localRadius;
	float volume;
	float viscosityConstant;
	float gasConstant;
	float restDensity;
	float particleMass;

	void instantiateFromFile(string file);

public:
	FluidSimulation(string file);
	void elapseTimeNaive();
	vector<Particle>& particleList();
};

#endif