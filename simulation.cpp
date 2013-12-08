#include <iostream>
#include <vector>

#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Particle
public:
	float massDensity;
	float pressure;
	vec3 position;
	vec3 velocity;
	vec3 acceleration;
};

class FluidSimulation {
private:
	vector<Particle> particles;

	vec3 gravity;
	float timeStepSize;
	int numGrids;
	int localRadius;
	float viscosityConstant;
	float gasConstant;
	float particleMass;
	float restDensity;
	float mass;

	int timeStep = 0;

	void instantiateFromFile(string file) {
		//store variables and set stuff at the end
		int width, height;

		ifstream inpfile(file.c_str());
		if(!inpfile.is_open()) {
			cout << "Unable to open file" << endl;
		} else {
			std::string line;
			//MatrixStack mst;

			while(inpfile.good()) {
				std::vector<std::string> splitline;
				std::string buf;

				std::getline(inpfile,line);
				std::stringstream ss(line);

				while (ss >> buf) {
					splitline.push_back(buf);
				}

				//Ignore blank lines
				if(splitline.size() == 0) {
					continue;
				}

				//Ignore comments
				if(splitline[0][0] == '#') {
					continue;
				}

				//Valid commands:
				//size width height
				//  must be first command of file, controls image size
				else if(!splitline[0].compare("size")) {
					width = atoi(splitline[1].c_str());
					height = atoi(splitline[2].c_str());
				}

				inpfile.close();
			}
		}
	}

public:
	FluidSimulation(string file) {

	}

	void elapseTimeNaive() {
		// Compute the density and pressure of each particle
		for (size_t i = 0; i < particles.size(); i++) {
			Particle& current = particles[i];
			current.massDensity = 0.0f;

			for (size_t j = 0; j < particles.size(); j++) {
				Particle& other = particles[j];
				current.density += other.mass * gaussianSmoothingKernel(current.position, other.position);
			}

			current.pressure = gasConstant * (pow(current.density / restDensity, 7.0f) - 1.0f);
		}

		// Update the position, velocity, and acceleration for each particle to the next time step
		for (size_t i = 0; i < particles.size(); i++) {
			Particle& current = particles[i];
			vec3 force;
			vec3 newAcceleration; // Acceleration at t+1

			// Advance position to time t+1
			// x_i+1 = x_i + v_i * delta_t + 0.5 * a_i * delta_t ^ 2
			current.position = current.position + current.velocity * timeStepSize + 0.5f * current.acceleration + pow(timeStepSize, 2.0f);

			// Force from gravity
			force = current.massDensity * gravity;

			for (size_t j = 0; j < particles.size(); j++) {
				Particle& other = particles[j];

				// Force from pressure
				force = force - force_pressure(current, other);

				// Force from viscosity
				force = force + viscosityConstant * force_viscosity(current, other);
			}

			newAcceleration = force / current.mass;

			// Advance velocity to time t+1
			// v_i+1 = v_i + 0.5 * (a_i + a_i+1) * delta_t
			current.velocity = current.velocity + 0.5f * (current.acceleration + newAcceleration) * timeStepSize;

			// Advance acceleration to time t+1
			current.acceleration = newAcceleration;	
		}
	}

	void writeToFile() {

	}

	vector<Particle>& particles() {
		return &particles;
	}
}