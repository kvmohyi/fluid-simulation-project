
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include "geometry.hpp"
#define PI 3.14159265f
using namespace std;
using namespace glm;

float sqr(float x) {
  return x * x;
}
/*
vec3 gradient_gaussian_smoothing(Particle particle, Particle& other){
	float r = glm::distance(particle.position, other.position);
	if (r > radius) {
		return vec3(0.0f, 0.0f, 0.0f);
	}

	vec3 vectorR  = vec3(other.position.x - particle.position.x, other.position.y - particle.position.y, other.position.z - other.position.z);
	vec3 gaussian = (pow(radius - r, 2) / r) * vectorR ;
	float constant = (-45 / (PI * pow(radius, 6)));
	gaussian = constant * gaussian;
	
	/*
	float gaussian = gaussian_smoothing(particle, other);  
	float partialX = -8 * abs(particle.position.x - other.position.x);
	float partialY = -8 * abs(particle.position.y - other.position.y);
	float partialZ = -8 * abs(particle.position.z - other.position.z);
	vec3 gradient = new vec3(partialX * gaussian, partialY * gaussian, partialZ * gaussian);*/
	/*return gaussian;

}

float gradient2_gaussian_smoothing(Particle particle, Particle& other){	
	float r = glm::distance(particle.position, other.position);
	if (r > radius) {
		return 0;
	}

	float gaussian = radius - r;
	gaussian = gaussian * 45 / (PI * pow(radius, 6));
	/*
	float gaussian = gaussian_smoothing(particle, other);
	float partialX = -8 * gaussian + 64 * sqr(abs(particle.position.x - other.position.x)) * gaussian;
	float partialY = -8 * gaussian + 64 * sqr(abs(particle.position.y - other.position.y)) * gaussian;
	float partialZ = -8 * gaussian + 64 * sqr(abs(particle.position.z - other.position.z)) * gaussian;
	vec3 gradient2 = new vec3(partialX, partialY, partialZ);
	*/
	/*return gaussian;
}
	
vec3 force_pressure(Particle particle, Particle& other, float mass) {
	float avgPressure = (particle.pressure + other.pressure) / 2.0;
	avgPressure = avgPressure * mass / other.massDensity;
	vec3 force = gradient_gaussian_smoothing(particle, other);
	force = avgPressure * force;
	return force;
}

vec3 force_viscosity(Particle particle, Particle& other, float mass) {
vec3 diffVelocity = other.velocity - particle.velocity;
diffVelocity = (mass / other.massDensity) * diffVelocity;
vec3 force = gradient2_gaussian_smoothing(particle, other) * diffVelocity;
return force;
}*/

float particleDistance(Particle& current, Particle& other) {
	return sqrt(sqr(current.position.x - other.position.x)
			  + sqr(current.position.y - other.position.y)
			  + sqr(current.position.z - other.position.z));
}

float spikyKernel(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	return 15.0f / (PI * pow(h, 6.0f)) * pow((h - r), 3.0f);
}

// Use this kernel for pressure
vec3 spikyKernelGradient(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	// it is current.position - other.position according to kelager p.18
	return -45.0f * (current.position - other.position) / (PI * pow(h, 6.0f) * r) * pow(h - r, 2.0f);
}

float viscosityKernel(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	return 15.0f / (2.0f * PI * pow(h, 3.0f)) * (pow(r, 3.0f) / (-2.0f * pow(h, 3.0f)) + pow(r, 2.0f) / pow(h, 2.0f) + h / (2.0f * r) - 1.0f);
}

// Use this kernel for viscosity
float viscosityKernelLaplacian(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	return 45.0f / (PI * pow(h, 6.0f)) * (h - r);
}

// The default kernel
float poly6Kernel(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	//cout << 315.0f / (64.0f * PI * pow(h, 9.0f)) * pow(pow(h, 2.0f) - pow(r, 2.0f), 3.0f) << endl;

	return 315.0f / (64.0f * PI * pow(h, 9.0f)) * pow(pow(h, 2.0f) - pow(r, 2.0f), 3.0f);
}

vec3 poly6KernelGradient(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	return -945.0f / (32.0f * PI * pow(h, 9.0f)) * (current.position - other.position) * pow((h * h - r * r), 2.0f);
}

float poly6KernelLaplacian(Particle& current, Particle& other, float h) {
	float r = particleDistance(current, other);

	return -945.0f / (32.0f * PI * pow(h, 9.0f)) * (h * h - r * r) * (3.0f * h * h - 7.0f * r * r);
}

vec3 pressureForcePartial(Particle& current, Particle& other, float particleMass, float h) {
	return -1.0f * particleMass * (current.pressure + other.pressure) / (2.0f * other.density) * spikyKernelGradient(current, other, h);
}

vec3 viscosityForcePartial(Particle& current, Particle& other, float timeStepSize, float viscosityConstant, float particleMass, float h) {

	return viscosityConstant * particleMass * (other.velocity - current.velocity) / other.density * viscosityKernelLaplacian(current, other, h);
}

vec3 inwardSurfaceNormalPartial(Particle& current, Particle& other, float particleMass, float h) {
	return particleMass / other.density * poly6KernelGradient(current, other, h);
}

float colorFieldLaplacianPartial(Particle& current, Particle& other, float particleMass, float h) {
	return particleMass / other.density * poly6KernelLaplacian(current, other, h);
}