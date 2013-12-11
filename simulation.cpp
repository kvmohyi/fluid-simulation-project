#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility> 

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <omp.h>

#include "simulation.hpp"
#include "equations.hpp"
#include "geometry.hpp"

using namespace std;
using namespace glm;

float dampFactor = 0.5;

void FluidSimulation::instantiateFromFile(string file) {
	std::ifstream inpfile(file.c_str());
	if(!inpfile.is_open()) {
		std::cout << "Unable to open file" << std::endl;
	} else {
		std::string line;

		while(inpfile.good()) {
			std::vector<std::string> splitline;
			std::string buf;

			std::getline(inpfile,line);
			std::stringstream ss(line);

			while (ss >> buf) {
				splitline.push_back(buf);
			}

			if(splitline.size() == 0) {
				continue;
			}
			else if(splitline[0][0] == '#') {
				continue;
			}
			else if(!splitline[0].compare("step_size")) {
				timeStepSize = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("num_grids")) {
				numGrids = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("grid_size")) {
				gridSize = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("radius")) {
				localRadius = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("world_size")) {
				worldSize = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("num_particles")) {
				numParticles = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("rest_density")) {
				restDensity = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("viscosity_constant")) {
				viscosityConstant = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("gas_constant")) {
				gasConstant = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("test_version")) {
				testVersion = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("dimensions")) {
				dimensions = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("tension_threshold")) {
				tensionThreshold = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("tension_constant")) {
				tensionConstant = atof(splitline[1].c_str());
			}
		}
		inpfile.close();
	}
}

FluidSimulation::FluidSimulation(string file) {
	instantiateFromFile(file);

	// Number of grids such that each grid cell is greater than the 2 * radius of support
	numGrids = ceil(worldSize / (2.0 * localRadius));
	// gridSize > 2 * localRadius
	gridSize = worldSize / numGrids;
	// Instantiate the "3D" grid of cells
	gridCells.resize(numGrids * numGrids * numGrids);
	// Set gravity here
	// Initialize the iteration count
	numIterations = 0;
	numParticles = 0;
	// Set up the test case
	drawTest(dimensions, testVersion);
	cube = RigidBody(worldSize * .75, worldSize * .75, worldSize * .75);
	//cout << cube.length << endl;
	// Verify the parameters
	printParams();
}

void FluidSimulation::elapseTimeGrid() {
	vector<vector<Particle> > newGridCells(numGrids * numGrids * numGrids);
	// Compute the density and pressure at each particle
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			current.density = 0.0f;

			int offsets[] = {-1, 0, 1};

			for (int z_offset = 0; z_offset < 3; z_offset++) {
				for (int y_offset = 0; y_offset < 3; y_offset++) {
					for (int x_offset = 0; x_offset < 3; x_offset++) {
						int index = mapToIndex(current, offsets[x_offset], offsets[y_offset], offsets[z_offset]);
						// If the cell is out of bounds, then ignore it
						if (index < 0 || index >= gridCells.size())
							continue;

						for (size_t j = 0; j < gridCells[index].size(); j++) {
							Particle& other = gridCells[index][j];

							float r = particleDistance(current, other);

							if (r > localRadius)
								continue;

							cout << "DERP" << endl;
							current.density += particleMass * poly6Kernel(current, other, localRadius);
						}
					}
				}
			}

			current.pressure = gasConstant * (current.density - restDensity);
			#if false
				cout << current.density << endl;
			#endif
		}
	}

	// Update the position, velocity, and acceleration for each particle to the next time step
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			vec3 pressureForce(0.0f, 0.0f, 0.0f);
			vec3 viscosityForce(0.0f, 0.0f, 0.0f);
			vec3 surfaceTensionForce(0.0f, 0.0f, 0.0f);
			vec3 inwardSurfaceNormal(0.0f, 0.0f, 0.0f);
			float colorFieldLaplacian = 0.0f;
			vec3 acceleration(0.0f, 0.0f, 0.0f);
			bool collide = false;
			float time = timeStepSize;

			#if true
				current.print();
				cout << endl;
			#endif

			int offsets[] = {-1, 0, 1};

			for (int z_offset = 0; z_offset < 3; z_offset++) {
				for (int y_offset = 0; y_offset < 3; y_offset++) {
					for (int x_offset = 0; x_offset < 3; x_offset++) {
						int index = mapToIndex(current, offsets[x_offset], offsets[y_offset], offsets[z_offset]);
						// If the cell is out of bounds, then ignore it
						if (index < 0 || index >= gridCells.size())
							continue;

						//cout << index << endl;
						//cout << "segfault" << endl;

						for (size_t j = 0; j < gridCells[index].size(); j++) {
							Particle& other = gridCells[index][j];

							// if i == j, then skip it
							if (index == i_cell && i == j)
								continue;
							
							float r = particleDistance(current, other);

							if (r > localRadius)
								continue;

							// Force from pressure
							pressureForce += pressureForcePartial(current, other, particleMass, localRadius);
							//cout << "Pressure Force: " << pressureForce.x << ", " << pressureForce.y << ", " << pressureForce.z << endl;
							// Force from viscosity
							viscosityForce += viscosityForcePartial(current, other, timeStepSize, viscosityConstant, particleMass, localRadius);
							// Surface tension yay!
							inwardSurfaceNormal += inwardSurfaceNormalPartial(current, other, particleMass, localRadius);
							colorFieldLaplacian += colorFieldLaplacianPartial(current, other, particleMass, localRadius);
						}
					}
				}
			}

			float inwardSurfaceNormalMagnitude = glm::length(inwardSurfaceNormal);

			cout << "Force: " << glm::to_string(pressureForce) << endl;

			if (inwardSurfaceNormalMagnitude >= tensionThreshold) {
				surfaceTensionForce = -1.0f * tensionConstant * colorFieldLaplacian * inwardSurfaceNormal / inwardSurfaceNormalMagnitude;
				//cout << "Inward Surface Tension Force: " << to_string(surfaceTensionForce) << endl;
			}
			// else leave it as the zero vector

			acceleration = (current.density * gravity + pressureForce /*+ viscosityForce + surfaceTensionForce*/) / current.density; // a at t=0
			
			/*vec3 oldVelocity;
			if (numIterations > 0)
				oldVelocity = current.currentVelocity;
			else
				oldVelocity = -0.5f * timeStepSize * acceleration;

			current.prevVelocity = oldVelocity + timeStepSize * acceleration;

			current.position = current.position + timeStepSize * current.prevVelocity;*/

			vec3 nextPrevVelocity;
			vec3 nextCurrentVelocity;
			vec3 nextPosition;

			if (numIterations == 0)
				current.prevVelocity = -0.5f * timeStepSize * acceleration;

			nextPrevVelocity = current.prevVelocity + timeStepSize * acceleration;
			nextCurrentVelocity = 0.5f * (current.prevVelocity + nextPrevVelocity);
			nextPosition = current.position + timeStepSize * nextPrevVelocity;

			Particle newParticle(nextPosition, nextCurrentVelocity, nextPrevVelocity, current.density, current.pressure);

			if (cube.collision(current.position, newParticle.position)) {
			  //cout << "collision"  << endl;
			  pair<float, vec3> timeNormal = cube.collisionTimeNormal(current.position, newParticle.position);
			  float timeStep = timeStepSize * timeNormal.first;

			  vec3 collisionPoint = current.position + timeNormal.first * (newParticle.position - current.position);

			  Particle collisionParticle = Particle(collisionPoint);

			  float depth = sqrt(sqr(newParticle.position.x - collisionParticle.position.x)
			  + sqr(newParticle.position.y - collisionParticle.position.y)
			  + sqr(newParticle.position.z - collisionParticle.position.z));

			  //cout << "depth" << depth << endl;
			  float rayDistance = sqrt(sqr(newParticle.position.x - current.position.x)
			  + sqr(newParticle.position.y - current.position.y)
			  + sqr(newParticle.position.z - current.position.z));

			  //cout << "rayDistance" << rayDistance << endl;
			  vec3 newVelocity = current.currentVelocity - (1 + dampFactor * depth / (rayDistance)) * dot(current.currentVelocity, timeNormal.second) * timeNormal.second;
			  //cout << "Collision Velocity " << newVelocity.x << " " << newVelocity.y << " " << newVelocity.z << endl;;
			  newParticle.currentVelocity = newVelocity;
			  newParticle.position = collisionPoint + (timeStepSize - timeStep) * newParticle.currentVelocity;
			}

			int newIndex = mapToIndex(newParticle);
			if (newIndex < 0 || newIndex >= gridCells.size())
				continue;
			newGridCells[mapToIndex(current)].push_back(newParticle);
		}
	}

	gridCells = newGridCells;

	numIterations++;
}

vector<vector<Particle> >& FluidSimulation::particleList() {
	return gridCells;
}

int FluidSimulation::mapToIndex(Particle particle) {
	int x_bucket = (worldSize / 2.0f + particle.position.x) / worldSize * numGrids;
	int y_bucket = (worldSize / 2.0f + particle.position.y) / worldSize * numGrids;
	int z_bucket = (worldSize / 2.0f + particle.position.z) / worldSize * numGrids;
	return mapToIndex(x_bucket, y_bucket, z_bucket);
}

int FluidSimulation::mapToIndex(Particle particle, int x_offset, int y_offset, int z_offset) {
	int x_bucket = x_offset + (worldSize / 2.0f + particle.position.x) / worldSize * numGrids;
	int y_bucket = y_offset + (worldSize / 2.0f + particle.position.y) / worldSize * numGrids;
	int z_bucket = z_offset + (worldSize / 2.0f + particle.position.z) / worldSize * numGrids;
	return mapToIndex(x_bucket, y_bucket, z_bucket);
}

int FluidSimulation::mapToIndex(int x, int y, int z) {
	return x + numGrids * (y + numGrids * z);
}

void FluidSimulation::drawWaterShape(int numParticles, float xStart, float yStart, float zStart, float xEnd, float yEnd, float zEnd){
	float length = xEnd - xStart;
	float width = yEnd- yStart;
	float depth = zEnd - zStart;

	if(depth == 0.0){
		float area = length * width;
		particleMass = area * restDensity / numParticles;
		float spacing = sqrt(numParticles / area);
		float step = 1 / spacing;		
		for(float x = xStart; x < xEnd; x += step)
		{
			for(float y = yStart; y < yEnd; y += step)
			{
				vec3 position = vec3(x, y, 0);
				Particle particle = Particle(position);
				int index = mapToIndex(particle);
				gridCells[index].push_back(particle);
			}
		}
	}						
	else {		
		float volume = length * width * depth;
		particleMass = volume * restDensity / numParticles;
		float spacing = pow(numParticles / volume, 1.0f / 3.0f);
		float step = 1/spacing;
		for(float x = xStart; x < xEnd; x += step)
		{	
			for(float y = yStart; y < yEnd; y += step)
			{	
				for(float z = zStart; z < zEnd; z += step)
				{
					vec3 position = vec3(x, y, z);
					Particle particle = Particle(position);
					int index = mapToIndex(particle);
					gridCells[index].push_back(particle);
				}
			}
		}
	}

	numIterations++;
}
void FluidSimulation::drawTest(int dimension, int version) {
	gravity = vec3(0.0f, -20.0f, 0.0);
	if(version == 1){
	  	particleMass = 0.2f;
		numParticles = 1;

		vec3 position1(0.0f, 0.0f, 0.0f);

		Particle particle1(position1);
		gridCells[mapToIndex(particle1)].push_back(particle1);
	}
	else if (version == 2) { // two particles
		particleMass = 0.2f;
		numParticles = 2;

		vec3 position1(-0.0130697f, 0.0f, 0.0f);
		vec3 position2(0.0130697f, 0.0f, 0.0f);

		Particle particle1(position1);
		Particle particle2(position2);

		gridCells[mapToIndex(particle1)].push_back(particle1);
		gridCells[mapToIndex(particle2)].push_back(particle2);
	}
	else if (version == 3) {
		particleMass = 0.02f;
		float stepSize = 0.025;
		for (float x = sideMax / -16.0f; x < sideMax / 16.0f; x += stepSize) {
			for (float y = sideMax / -16.0f; y < sideMax / 16.0f; y += stepSize) {
				for (float z = sideMax / -16.0f; z < sideMax / 16.0f; z += stepSize) {
					vec3 position(x, y, z);
					Particle particle(position);
					gridCells[mapToIndex(particle)].push_back(particle);
					numParticles++;
				}
			}
		}
	}
}

void FluidSimulation::printParams() {
	cout << "gravity: <" << gravity.x << ", " << gravity.y << ", " << gravity.z << ">" << endl;
	cout << "timeStepSize: " << timeStepSize << endl;
	cout << "numGrids: " << numGrids << endl;
	cout << "gridCells.size(): " << gridCells.size() << endl;
	cout << "numPartices: " << numParticles << endl;
	cout << "worldSize: " << worldSize << endl;
	cout << "gridSize: " << gridSize << endl;
	cout << "localRadius: " << localRadius << endl;
	cout << "viscosityConstant: " << viscosityConstant << endl;
	cout << "gasConstant: " << gasConstant << endl;
	cout << "restDensity: " << restDensity << endl;
	cout << "particleMass: " << particleMass << endl;
	cout << "testVersion: " << testVersion << endl;
	cout << "dimensions: " << dimensions << endl;
	cout << "tensionThreshold: " << tensionThreshold << endl;
	cout << "tensionConstant: " << tensionConstant << endl;
}