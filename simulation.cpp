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

using namespace std;
using namespace glm;

float dampFactor = 0.5f;
#define PI 3.14159265f

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
			else if(!splitline[0].compare("kernel_x")) {
				kernelX = atof(splitline[1].c_str());
			}
		}
		inpfile.close();
	}
}

FluidSimulation::FluidSimulation(string file) {
	gravity = vec3(0.0f, -9.8f, 0.0);
	worldSize = 2.0f;
	timeStepSize = 0.01;
	numIterations = 0;
	//numParticles = 0;
	//localRadius = 0.0625;
	restDensity = 998.29;
	viscosityConstant = 3.5;
	gasConstant = 3.0;
	tensionConstant = 0.0728;
	tensionThreshold = 7.065;
	testVersion = 5;
	dimensions = 3;
	//Set num grids from the test cases
	//numGrids = ceil(worldSize / (2.0f * localRadius));
	//gridSize = worldSize / numGrids;
	//gridCells.resize(numGrids * numGrids * numGrids);
	instantiateFromFile(file);
	drawTest(dimensions, testVersion);
	cube = RigidBody(worldSize, worldSize, worldSize);
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

			#if false
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
							
							float r = particleDistance(current, other);

							if (r > localRadius)
								continue;

							// Surface tension yay!
							inwardSurfaceNormal += inwardSurfaceNormalPartial(current, other, particleMass, localRadius);
							colorFieldLaplacian += colorFieldLaplacianPartial(current, other, particleMass, localRadius);

							// if i == j, then skip it
							if (index == i_cell && i == j)
								continue;

							// Force from pressure
							pressureForce += pressureForcePartial(current, other, particleMass, localRadius);
							//cout << "Pressure Force: " << pressureForce.x << ", " << pressureForce.y << ", " << pressureForce.z << endl;
							// Force from viscosity
							viscosityForce += viscosityForcePartial(current, other, timeStepSize, viscosityConstant, particleMass, localRadius);
						}
					}
				}
			}

			float inwardSurfaceNormalMagnitude = glm::length(inwardSurfaceNormal);

			//cout << "\tForce: " << glm::to_string(pressureForce) << endl;

			if (inwardSurfaceNormalMagnitude >= tensionThreshold) {
				surfaceTensionForce = -1.0f * tensionConstant * colorFieldLaplacian * inwardSurfaceNormal / inwardSurfaceNormalMagnitude;
				//cout << "Inward Surface Tension Force: " << to_string(surfaceTensionForce) << endl;
			}
			// else leave it as the zero vector

			acceleration = (current.density * gravity + pressureForce /*+ viscosityForce + surfaceTensionForce*/) / current.density; // a at t=0
			
			/*vec3 oldVelocity;
			if (numIterations > 0)
				oldVelocity = current.velocity;
			else
				oldVelocity = -0.5f * timeStepSize * acceleration;

			current.prevVelocity = oldVelocity + timeStepSize * acceleration;

			current.position = current.position + timeStepSize * current.prevVelocity;*/

			/*vec3 nextPrevVelocity;
			vec3 nextCurrentVelocity;
			vec3 nextPosition;

			if (numIterations == 0)
				current.prevVelocity = -0.5f * timeStepSize * acceleration;

			nextPrevVelocity = current.prevVelocity + timeStepSize * acceleration;
			nextCurrentVelocity = 0.5f * (current.prevVelocity + nextPrevVelocity);
			nextPosition = current.position + timeStepSize * nextPrevVelocity;*/

			vec3 nextPrevVelocity(0.0f, 0.0f, 0.0f);
			vec3 nextCurrentVelocity = current.velocity + timeStepSize * acceleration;
			vec3 nextPosition = current.position + timeStepSize * (current.velocity + timeStepSize * acceleration);

			Particle newParticle(nextPosition, nextCurrentVelocity, current.density, current.pressure);

			/*if (cube.collision(current.position, newParticle.position)) {
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
			  vec3 newVelocity = current.velocity - (1 + dampFactor * depth / (rayDistance)) * dot(current.velocity, timeNormal.second) * timeNormal.second;
			  //cout << "Collision Velocity " << newVelocity.x << " " << newVelocity.y << " " << newVelocity.z << endl;;
			  newParticle.velocity = newVelocity;
			  newParticle.position = collisionPoint + (timeStepSize - timeStep) * newParticle.velocity;
			  //newParticle.prevVelocity = newParticle.velocity - 0.5f * timeStep * acceleration;
			  //cout << "newPosition " << newParticle.position.x << " " << newParticle.position.y << " " << newParticle.position.z << endl;
			}*/

			 cube.handleCollisions(newParticle);

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
  

void FluidSimulation::drawWaterShape(int n, float xStart, float yStart, float zStart, float xEnd, float yEnd, float zEnd, vec3 velocity){
	float length = xEnd - xStart;
	float width = yEnd- yStart;
	float depth = zEnd - zStart;

	if(abs(depth) < 0.1){
		float area = length * width;
		//particleMass = area * restDensity / n;
		float ratio = length / width;
		float particlesWidth = pow((float)n / ratio, 1.0f / 2.0f);
		float spacing = width / particlesWidth;
		for(float x = xStart; x < xEnd; x += spacing)
		{
			for(float y = yStart; y < yEnd; y += spacing)
			{
				vec3 position = vec3(x, y, 0);
				Particle particle = Particle(position, velocity);
				int index = mapToIndex(particle);
				gridCells[index].push_back(particle);
			}
		}
	}						
	else {		
		float volume = length * width * depth;
		//particleMass = volume * restDensity / n;
		float ratioY = width / length;
		float ratioZ = depth / length;
		float particlesLength = pow((float) n / (ratioY * ratioZ), 1.0f / 3.0f);
		float spacing = length / particlesLength;
		for(float x = xStart; x < xEnd; x += spacing)
		{	
			for(float y = yStart; y < yEnd; y += spacing)
			{	
				for(float z = zStart; z < zEnd; z += spacing)
				{
					vec3 position = vec3(x, y, z);
					Particle particle = Particle(position, velocity);
					int index = mapToIndex(particle);
					gridCells[index].push_back(particle);
				}
			}
		}
	}

	numIterations++;
}
void FluidSimulation::drawTest(int dimension, int version) {
	float sideMax = worldSize / 2.0f;

	if (version == 0) { // One particle
	  	particleMass = 0.2f;
		numParticles = 1;

		vec3 position1(0.0f, 0.0f, 0.0f);

		Particle particle1(position1);
		gridCells[mapToIndex(particle1)].push_back(particle1);
	}
	else if (version == 1){
	  vec3 initVelocity = vec3(0.0, 0.0, 0.0);
	  if (dimensions == 3) {
	    float volume = pow(sideMax, 3.0);
	    particleMass = restDensity * volume / numParticles;
	    localRadius = idealLocalRadius3D(volume);
	    reinitGridCells();
	    drawWaterShape(numParticles, -1.0 * sideMax, -1.0 * sideMax, -1.0 * sideMax, 0.0, 0.0, 0.0, initVelocity);
	  }
	  else {
	    float area = pow(sideMax, 2.0);
	    particleMass = restDensity * area / numParticles;
	    localRadius = idealLocalRadius2D(area);
	    reinitGridCells();
	    drawWaterShape(numParticles, -1.0 * sideMax, -1.0 * sideMax, 0.0, 0.0, 0.0, 0.0, initVelocity);
	  }
	}
	else if (version == 2) {
		if (dimensions == 3) {
			float cubeSize = worldSize / 3.0;
			float poolXZSize = worldSize;
			float poolYSize = worldSize / 8.0;
			vec3 initVelocity = vec3(0.0f, 0.0f, 0.0f);

			float cubeVolume = cubeSize * cubeSize * cubeSize;
			float poolVolume = poolXZSize * poolXZSize * poolYSize;
			float totalVolume = cubeVolume + poolVolume;
			int particlesInCube = cubeVolume / totalVolume * numParticles;
			int particlesInPool = poolVolume / totalVolume * numParticles;

			localRadius = idealLocalRadius3D(totalVolume);
			particleMass = totalVolume * restDensity / numParticles;
			reinitGridCellls();

			// Draw the cube
			drawWaterShape(particlesInCube, -cubeSize / 2.0, worldSize / 6.0, -cubeSize / 2.0, cubeSize / 2.0, worldSize / 2.0, cubeSize / 2.0, initVelocity);
			// Draw the pool
			drawWaterShape(particlesInPool, -worldSize / 2.0, -worldSize / 2.0, -worldSize / 2.0, worldSize / 2.0, -worldSize / 2.0 + poolYSize / 2.0, worldSize / 2.0, initVelocity);
		}
		else { // dimensions == 2
			float cubeSize = worldSize / 3.0;
			float poolXSize = worldSize;
			float poolYSize = worldSize / 8.0;
			vec3 initVelocity = vec3(0.0f, 0.0f, 0.0f);

			float cubeArea = cubeSize * cubeSize;
			float poolArea = poolXSize * poolYSize;
			float totalArea = cubeArea + poolArea;
			int particlesInCube = cubeArea / totalArea * numParticles;
			int particlesInPool = poolArea / totalArea * numParticles;

			localRadius = idealLocalRadius2D(totalArea);
			particleMass = totalArea * restDensity / numParticles;
			reinitGridCellls();

						// Draw the cube
			drawWaterShape(particlesInCube, -cubeSize / 2.0, worldSize / 6.0, 0.0, cubeSize / 2.0, worldSize / 2.0, 0.0, initVelocity);
			// Draw the pool
			drawWaterShape(particlesInPool, -worldSize / 2.0, -worldSize / 2.0, 0.0, worldSize / 2.0, -worldSize / 2.0 + poolYSize / 2.0, 0.0, initVelocity);
		}
	}
	else if(version == 3) {
	  vec3 initVelocity1 = vec3(-1.0f, 0.0f, 0.0f);
	  vec3 initVelocity2 = vec3(1.0f, 0.0f, 0.0f);
	  if(dimensions == 3){
	    float volume = pow(worldSize / 3.0f, 3.0) * 2;
	    particleMass = restDensity * volume / numParticles;
	    localRadius = idealLocalRadius3D(volume);
	    reinitGridCells();
	    drawWaterShape(numParticles / 2.0f, -1.0 * sideMax, 0.5 * sideMax, -0.25 * sideMax, -0.5 * sideMax, 1.0 * sideMax, 0.25 * sideMax, initVelocity2);
	    drawWaterShape(numParticles / 2.0f, 0.5 * sideMax, 0.5 * sideMax, -0.25 * sideMax, 1.0 * sideMax, 1.0 * sideMax, 0.25 * sideMax, initVelocity1);
	  }
	  else {
	    float area = pow(worldSize / 3.0f, 2.0) * 2;
	    particleMass = restDensity * area / numParticles;
	    localRadius = ideaLocalRadius2D(area);
	    reinitGridCells();
	    drawWaterShape(numParticles / 2.0f, -1.0 * sideMax, 0.5 * sideMax, 0.0, -0.5 * sideMax, 1.0 * sideMax, 0.0, initVelocity2);
	    drawWaterShape(numParticles / 2.0f, 0.5 * sideMax, 0.5 * sideMax, 0.0, 1.0 * sideMax, 1.0 * sideMax, 0.0, initVelocity1);
	  }
	  	    
	}
	else if (version == 6) {
		float sideSize = 1.0;
		float volume = sideSize * sideSize * sideSize;
		float maxParticles = 100000;
		float stepSize = sideSize / pow(maxParticles, 1.0 / 3.0);
		particleMass = volume * restDensity / maxParticles;
		localRadius = idealLocalRadius3D(volume);
		reinitGridCellls();
		for (float x = -worldSize / 2.0; x < 0.0; x += stepSize) {
			for (float y = -worldSize / 2.0; y < 0.0; y += stepSize) {
				for (float z = -worldSize / 2.0; z < 0.0; z += stepSize) {
					vec3 position(x, y, z);
					Particle particle(position);
					gridCells[mapToIndex(particle)].push_back(particle);
					numParticles++;
				}
			}
		}
	}
	/*else if (version == 3) {
		float volume = 0.1f;
		float sideSize = pow(volume, 1.0f / 3.0f);

		restDensity = 1000;
		numParticles = 5000;
		particleMass = volume * restDensity / numParticles;

		//the smallest possible amount of particles that renders the fluid simulation stable,
		//while still respecting the properties of the fluid material.
		float x = 20;
		//localRadius = pow((3.0f * volume * x) / (4 * PI * numParticles), 1.0f / 3.0f);
		localRadius = idealLocalRadius(volume);
		float stepSize = pow(volume / numParticles, 1.0f / 3.0f);
		numParticles = 0;

		for (float x = sideSize / -2.0f + stepSize / 2.0f; x < sideSize / 2.0f; x += stepSize) {
			for (float y = sideSize / -2.0f + stepSize / 2.0f; y < sideSize / 2.0f; y += stepSize) {
				for (float z = sideSize / -2.0f + stepSize / 2.0f; z < sideSize / 2.0f; z += stepSize) {
					vec3 position(x, y, z);
					Particle particle(position);
					gridCells[mapToIndex(particle)].push_back(particle);
					numParticles++;
				}
			}
		}
	}
	else if (version == 4) { // cube drop
		float sideSize = 0.5f;

		particleMass = 0.02;
		localRadius = 0.05;
		reinitGridCells();
		float stepSize = 0.01;

		for (float x = sideSize / -2.0f; x < sideSize / 2.0f; x += stepSize) {
			for (float y = sideSize / -2.0f; y < sideSize / 2.0f; y += stepSize) {
				vec3 position(x, y, 0.0f);
				Particle particle(position);
				gridCells[mapToIndex(particle)].push_back(particle);
				numParticles++;
				for (float z = sideSize / -2.0f; z < sideSize / 2.0f; z += stepSize) {
					
				}
			}
		}
	}*/
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
	cout << "kernelX" << kernelX << endl;
}

double FluidSimulation::sphereRadius(Particle& particle) {
	return pow(3.0 * particleMass / (4.0 * PI * particle.density), 1.0 / 3.0);
}

float FluidSimulation::idealLocalRadius3D(float volume) {
	return pow(3.0 * volume * kernelX / (4.0 * PI * numParticles), 1.0 / 3.0);
}

float FluidSimulation::idealLocalRadius2D(float area) {
	return pow(area * kernelX / (PI * numParticles), 1.0 / 2.0);
}

void FluidSimulation::reinitGridCells() {
	numGrids = ceil(worldSize / (2.0f * localRadius));
	gridSize = worldSize / numGrids;
	gridCells.resize(numGrids * numGrids * numGrids);
}
