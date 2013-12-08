#ifndef	_SIMULATION_HPP
#define _SIMULATION_HPP

#include <vector>

#include <glm/glm.hpp>

#include "geometry.hpp"

using namespace std;
using namespace glm;

class FluidSimulation {
private:
	vector<Particle> particles;
	vector<RigidBody> rigidBodies;

	vec3 gravity;
	float timeStepSize;
	int numGrids;
	int localRadius;
	float viscosityConstant;
	float gasConstant;
	float particleMass;
	float restDensity;

	void instantiateFromFile(string file);

public:
	FluidSimulation(string file);
	void elapseTimeNaive();
	vector<Particle>& particleList();
};

#endif