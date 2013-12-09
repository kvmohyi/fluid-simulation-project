
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <glm/glm.hpp>

#include "geometry.hpp"
#define PI 3.14159265
using namespace std;
using namespace glm;

/*float radius = 10.0;

float sqr(float x) {
  return x * x;
}
	    
float gaussian_smoothing(Particle particle, Particle other){
	float r = glm::distance(particle.position, other.position);
	if (r > radius) {
		return 0;
	}
	float gaussian = pow(sqr(radius) - sqr(r), 3);
	float constant = 315 / (64 * PI * pow(radius, 9)) ;
	gaussian = constant * gaussian;
	//	float gaussian = exp(-4 * sqr(r/(float)(radius));
	return gaussian;
}

vec3 gradient_gaussian_smoothing(Particle particle, Particle other){
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

float gradient2_gaussian_smoothing(Particle particle, Particle other){	
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
	
vec3 force_pressure(Particle particle, Particle other, float mass) {
	float avgPressure = (particle.pressure + other.pressure) / 2.0;
	avgPressure = avgPressure * mass / other.massDensity;
	vec3 force = gradient_gaussian_smoothing(particle, other);
	force = avgPressure * force;
	return force;
}

vec3 force_viscosity(Particle particle, Particle other, float mass) {
vec3 diffVelocity = other.velocity - particle.velocity;
diffVelocity = (mass / other.massDensity) * diffVelocity;
vec3 force = gradient2_gaussian_smoothing(particle, other) * diffVelocity;
return force;
}*/

// This does not include the mass multiplier
vec3 pressureForcePartial(Particle current, Particle other, float h) {
	return -1.0f * (current.density + other.density) / (2 * other.density) * gradient_gaussian_smoothing(current, other);
}

// This does not include the mass and viscosity constant multipliers
vec3 viscosityForcePartial(Particle current, Particle other, float h) {
	return (other.velocity - current.velocity) / other.density * gradient2_gaussian_smoothing(current, other);
}

// Use this one for presure
float spikyKernel(Particle current, Particle other, float h) {
	vec3 r = glm::distance(current, other);

	if (r < 0.0f || r > h)
		return 0.0f;

	return 15.0f / (PI * pow(h, 6.0f)) * pow((h - r), 3.0f);
}

float viscosityKernel(Particle current, Particle other, float h) {
	vec3 r = glm::distance(current, other);

	if (r < 0.0f || r > h)
		return 0.0f;

	return 15.0f / (2.0f * PI * pow(h, 3.0f)) * (pow(r, 3.0f) / (-2.0f * pow(h, 3.0f)) + pow(r, 2.0f) / pow(h, 2.0f) + h / (2.0f * r) - 1.0f);
}

float poly6Kernel(Particle current, Particle other, float h) {
	vec3 r = glm::distance(current, other);

	if (r < 0.0f || r > h)
		return 0.0f;

	return 315.0f / (64.0f * PI * pow(h, 9.0f)) * pow(pow(h, 2.0f) - pow(r, 2.0f), 3.0f);
}