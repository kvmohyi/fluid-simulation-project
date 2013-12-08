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
			else if(!splitline[0].compare("particle")) {

				// jonathan pls fix

				/*Particle particle = vec3(atof(splitline[1].c_str()), atof(splitline[2].c_str()), atof(splitline[3].c_str()));
				particles.push_back(particle);*/
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
			else if(!splitline[0].compare("radius")) {
				localRadius = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("num_particles")) {
				numParticles = atoi(splitline[1].c_str());
			}
			else if(!splitline[0].compare("volume")) {
				volume = atof(splitline[1].c_str());
			}
			else if(!splitline[0].compare("mass")) {
				particleMass = atof(splitline[1].c_str());
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
			else if(!splitline[0].compare("gravity")) {
				gravity = vec3(atof(splitline[1].c_str()), atof(splitline[2].c_str()), atof(splitline[3].c_str()));
			}
		}
		inpfile.close();
	}
}

FluidSimulation::FluidSimulation(string file) {
	instantiateFromFile(file);
}

void FluidSimulation::elapseTimeNaive() {
	// Compute the density and pressure of each particle
	for (size_t i = 0; i < particles.size(); i++) {
		Particle& current = particles[i];
		current.massDensity = 0.0f;

		for (size_t j = 0; j < particles.size(); j++) {
			Particle& other = particles[j];
			current.massDensity += particleMass * gaussian_smoothing(current, other);
		}

		current.pressure = gasConstant * (pow(current.massDensity * particleMass / restDensity, 7.0f) - 1.0f);
	}

	// Update the position, velocity, and acceleration for each particle to the next time step
	// using leapfrog integration
	for (size_t i = 0; i < particles.size(); i++) {
		Particle& current = particles[i];
		vec3 force;
		vec3 newPosition; // Position at t+1
		vec3 newAcceleration; // Acceleration at t+1

		// Advance position to time t+1
		// x_i+1 = x_i + v_i * delta_t + 0.5 * a_i * delta_t ^ 2
		newPosition = current.position + current.velocity * timeStepSize + 0.5f * current.acceleration + pow(timeStepSize, 2.0f);
		current.position = newPosition;

		// Force from gravity
		force = current.massDensity * gravity;

		for (size_t j = 0; j < particles.size(); j++) {
			Particle& other = particles[j];

			// Force from pressure
			force = force - force_pressure(current, other);

			// Force from viscosity
			force = force + viscosityConstant * force_viscosity(current, other);
		}

		newAcceleration = force / current.massDensity;

		// Advance velocity to time t+1
		// v_i+1 = v_i + 0.5 * (a_i + a_i+1) * delta_t
		current.velocity = current.velocity + 0.5f * (current.acceleration + newAcceleration) * timeStepSize;

		// Advance acceleration to time t+1
		current.acceleration = newAcceleration;	
	}
}

vector<Particle>& FluidSimulation::particleList() {
	return particles;
}