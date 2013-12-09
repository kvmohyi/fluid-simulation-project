#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include <glm/glm.hpp>

#include "simulation.hpp"
#include "equations.hpp"
#include "geometry.hpp"

using namespace std;
using namespace glm;

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
				gridSize = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("radius")) {
				localRadius = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("world_size")) {
				worldSize = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("num_particles")) {
				numParticles = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("rest_density")) {
				restDensity = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("viscosity_constant")) {
				viscosityConstant = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("gas_constant")) {
				gasConstant = atoi(splitline[1].c_str());
			}
		}
		inpfile.close();
	}
}

FluidSimulation::FluidSimulation(string file) {
	instantiateFromFile(file);
}

void FluidSimulation::elapseTimeGrid() {
	// Compute the density and pressure of each particle
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			current.massDensity = 0.0f;

			for (size_t  j_cell = 0; j_cell < gridCells.size(); j_cell++) {
				for (size_t j = 0; j < gridCells[j_cell].size(); j++) {
					Particle& other = gridCells[j_cell][j];
					current.massDensity += particleMass * gaussian_smoothing(current, other);
				}
			}

			current.pressure = gasConstant * (pow(current.massDensity * particleMass / restDensity, 7.0f) - 1.0f);
		}
	}

	// Update the position, velocity, and acceleration for each particle to the next time step
	// using leapfrog integration
	for (size_t  i_cell = 0; i_cell < gridCells.size(); i_cell++) {
		for (size_t i = 0; i < gridCells[i_cell].size(); i++) {
			Particle& current = gridCells[i_cell][i];
			vec3 force;
			vec3 newPosition; // Position at t+1
			vec3 newAcceleration; // Acceleration at t+1

			// Advance position to time t+1
			// x_i+1 = x_i + v_i * delta_t + 0.5 * a_i * delta_t ^ 2
			newPosition = current.position + current.velocity * timeStepSize + 0.5f * current.acceleration + pow(timeStepSize, 2.0f);
			current.position = newPosition;

			// Force from gravity
			force = current.massDensity * gravity;

			int offsets[] = {-1, 1};

			for (int x_offset = 0; x_offset < 2; x_offset++) {
				for (int y_offset = 0; y_offset < 2; y_offset++) {
					for (int z_offset = 0; z_offset < 2; z_offset++) {
						int index = mapToIndex(offsets[x_offset], offsets[y_offset], offsets[z_offset]);
						// If the cell is out of bounds, then ignore it
						if (index < 0 || index >= gridCells.size())
							continue;

						for (size_t j = 0; j < gridCells[index].size(); j++) {
							Particle& other = gridCells[index][j];
							// Force from pressure
							force = force - force_pressure(current, other);
							// Force from viscosity
							force = force + viscosityConstant * force_viscosity(current, other);
						}
					}
				}
			}

			for (size_t  j_cell = 0; j_cell < gridCells.size(); j_cell++) {
				
			}

			newAcceleration = force / current.massDensity;

			// Advance velocity to time t+1
			// v_i+1 = v_i + 0.5 * (a_i + a_i+1) * delta_t
			current.velocity = current.velocity + 0.5f * (current.acceleration + newAcceleration) * timeStepSize;

			// Advance acceleration to time t+1
			current.acceleration = newAcceleration;
		}
	}
}

vector<Particle>& FluidSimulation::particleList() {
	vector<Particle> asdf;
	return asdf;
}

int FluidSimulation::mapToIndex(Particle particle) {
	int x_bucket = (worldSize / 2 + particle.position.x) / worldSize * numGrids;
	int y_bucket = (worldSize / 2 + particle.position.y) / worldSize * numGrids;
	int z_bucket = (worldSize / 2 + particle.position.z) / worldSize * numGrids;
	return mapToIndex(x_bucket, y_bucket, z_bucket);
}

int FluidSimulation::mapToIndex(Particle particle, int x_offset, int y_offset, int z_offset) {
	int x_bucket = x_offset + (worldSize / 2 + particle.position.x) / worldSize * numGrids;
	int y_bucket = y_offset + (worldSize / 2 + particle.position.y) / worldSize * numGrids;
	int z_bucket = z_offset + (worldSize / 2 + particle.position.z) / worldSize * numGrids;
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
		float spacing = pow(numParticles / volume, .333333333);
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
}
void FluidSimulation::drawTest(int dimension, int version){
	if(dimension == 2)
	{
		if(version == 1){
		  drawWaterShape(numParticles, -2 / (worldSize / 2), -2 / (worldSize / 2), 0, 2 / (worldSize / 2), 2 / (worldSize / 2), 0);//just one water cube in 2D
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
		  drawWaterShape(numParticles, -2 / (worldSize / 2), -2 / (worldSize / 2),  -2 / (worldSize / 2) , 2 / (worldSize / 2), 2 / (worldSize / 2), 2 / (worldSize / 2));
		}
		/*
		else {
			drawWaterShape(num_particles * .2, xStart, yStart, zStart, xEnd, yEnd, zEnd);
			drawWaterShape(num_particles * .8,-1, -1, -1, 1, 1, 1);
			}*/
	} 
	
}
