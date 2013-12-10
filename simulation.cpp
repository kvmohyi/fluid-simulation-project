#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <utility> 

#include <glm/glm.hpp>

#include "simulation.hpp"
#include "equations.hpp"
#include "geometry.hpp"

using namespace std;
using namespace glm;

float dampFactor = 0.7;

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
		}
		inpfile.close();
	}
}

FluidSimulation::FluidSimulation(string file) {
	instantiateFromFile(file);

	// Number of grids such that each grid cell is greater than the radius of support
	numGrids = ceil(worldSize / localRadius);
	// gridSize > localRadius
	gridSize = worldSize / numGrids;
	// Instantiate the 3D grid of cells
	gridCells.resize(numGrids * numGrids * numGrids);
	// Set gravity here
	gravity = vec3(0.0f, 0.0f, 0.0f);
	// Initialize the iteration count
	numIterations = 0;
	// Set up the test case
	drawTest(dimensions, testVersion);
	cube = RigidBody((float)(worldSize) * 0.7, (float)(worldSize) * 0.7, (float)(worldSize) * 0.7);
	cout << cube.length << endl;
	// Verify program input
	printParams();
}

void FluidSimulation::elapseTimeGrid() {
	vector<vector<Particle> > newGridCells(numGrids * numGrids * numGrids);
	// Compute the density and pressure at each particle
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			current.density = 0.0f;

			for (size_t  j_cell = 0; j_cell < gridCells.size(); j_cell++) {
				for (size_t j = 0; j < gridCells[j_cell].size(); j++) {
					Particle& other = gridCells[j_cell][j];
					current.density += particleMass * poly6Kernel(current, other, localRadius);
				}
			}
			//current.pressure = gasConstant * (pow(current.massDensity * particleMass / restDensity, 7.0f) - 1.0f);

			current.pressure = gasConstant * (current.density - restDensity);
			#if true
			cout << current.pressure << endl;
			#endif
		}
	}

	// Update the position, velocity, and acceleration for each particle to the next time step
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			vec3 pressureForce(0.0f, 0.0f, 0.0f);
			vec3 viscosityForce(0.0f, 0.0f, 0.0f);
			vec3 newPosition; // Position at t+1
			vec3 newAcceleration; // Acceleration at t+1
			vec3 newVelocity;
			bool collide = false;
			float time = timeStepSize;

			// Advance position to time t+1 using leapfrog integration
			// x_i+1 = x_i + v_i * delta_t + 0.5 * a_i * delta_t ^ 2
			//newPosition = current.position + current.velocity * timeStepSize + 0.5f * current.acceleration + pow(timeStepSize, 2.0f);

			if(cube.collision(current.position, newPosition)){
			  collide = true;
			  cout << "collision" << endl;
			  pair<float, vec3> timeAndNormal = cube.collisionTimeNormal(current.position, newPosition);
			  cout << "normal " << timeAndNormal.second.x << " " << timeAndNormal.second.y << " " << timeAndNormal.second.z << endl;
			  cout << "dot " << dot(current.velocity, timeAndNormal.second) << endl;
			    vec3 currVelocity = current.velocity;
			    vec3 reflectVelocity = currVelocity - 2 * dot(currVelocity, timeAndNormal.second) * timeAndNormal.second;
			    cout << "reflective Velocity: " << reflectVelocity.x << " " << reflectVelocity.y << " " << reflectVelocity.z << endl;
			    current.velocity = reflectVelocity;
			    newVelocity = reflectVelocity;
			    float newTime = timeStepSize - timeAndNormal.first * timeStepSize;
			    time = newTime;
			    vec3 collideLoc = current.position + (newPosition - current.position) * timeAndNormal.first;
			    newPosition = collideLoc + dampFactor * reflectVelocity * time;// + 0.5f * current.acceleration + pow(time, 2.0f);
			    cout << "newPosition " << newPosition.x << " " << newPosition.y << " " << newPosition.z << endl;
			    
			}
			current.position = newPosition;

			//current.position = newPosition;

			int offsets[] = {-1, 1};

			for (int x_offset = 0; x_offset < 2; x_offset++) {
				for (int y_offset = 0; y_offset < 2; y_offset++) {
					for (int z_offset = 0; z_offset < 2; z_offset++) {
						int index = mapToIndex(current, offsets[x_offset], offsets[y_offset], offsets[z_offset]);
						// If the cell is out of bounds, then ignore it
						if (index < 0 || index >= gridCells.size())
							continue;

						//cout << index << endl;
						//cout << "segfault" << endl;

						for (size_t j = 0; j < gridCells[index].size(); j++) {
							Particle& other = gridCells[index][j];
							// Force from pressure
							pressureForce += pressureForcePartial(current, other, particleMass, localRadius);
							// Force from viscosity
							viscosityForce += viscosityForcePartial(current, other, viscosityConstant, particleMass, localRadius);
						}
					}
				}
			}

			//cout << "got here" << endl;

			//cout << "Pressure Force: " << pressureForce.x << ", " << pressureForce.y << ", " << pressureForce.z << endl;

			newAcceleration = (current.density * gravity + pressureForce + viscosityForce) / current.density;
			current.acceleration = new Acceleration;
			// Advance velocity to time t+1 using leapfrog integration
			// v_i+1 = v_i + 0.5 * (a_i + a_i+1) * delta_t
			//current.velocity = current.velocity + 0.5f * (current.acceleration + newAcceleration) * timeStepSize;

			// Advance acceleration to time t+1
			//current.acceleration = newAcceleration;

			//cout << "Acceleration: " << current.acceleration.x << ", " << current.acceleration.y << ", " << current.acceleration.z << endl;

			int newIndex = mapToIndex(current);
			if (newIndex < 0 || newIndex >= gridCells.size())
				continue;
			newGridCells[mapToIndex(current)].push_back(current);
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
void FluidSimulation::drawTest(int dimension, int version){
	if(dimension == 2)
	{
		if(version == 1){
			drawWaterShape(numParticles, -.3, worldSize / 4.0, 0, .3, worldSize / 2.0, 0);
		  //drawWaterShape(numParticles, -.3, -2 / (worldSize / 2), 0, 2 / (worldSize / 2), 2 / (worldSize / 2), 0);//just one water cube in 2D
		}
		/*
		else 
		{
		        drawWaterShape(num_particles * .2, xStart, yStart, 0, xEnd, yEnd, 0);// dropping a cube into a body of water
			drawWaterShape(num_particles * .8, -1, -1, 0, 1, 1, 0);
			}*/
	}
	else
	{
		if(version == 1){
		  drawWaterShape(numParticles, worldSize / -8.0, worldSize / 4.0, worldSize / -8.0 , worldSize / 8.0, worldSize / 2.0, worldSize / 8.0);
		}
		/*
		else {
			drawWaterShape(num_particles * .2, xStart, yStart, zStart, xEnd, yEnd, zEnd);
			drawWaterShape(num_particles * .8,-1, -1, -1, 1, 1, 1);
			}*/
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
}