
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

float radius = 0.0625;

float sqr(float x) {
  return x * x;
}

float particle_distance(Particle particle, Particle other){
        vec3 pos1 = particle.position;
	vec3 pos2 = other.position;
	return sqrt(pow(pos2.x - pos1.x, 2) + pow(pos2.y - pos1.y, 2) + pow(pos2.z - pos1.z, 2));
}		    
float gaussian_smoothing(Particle particle, Particle other){
	float r = particle_distance(particle, other);
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
	float r = particle_distance(particle, other);
	if (r > radius) {
		return vec3(0.0, 0.0, 0.0);
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
	return gaussian;

}

float gradient2_gaussian_smoothing(Particle particle, Particle other){	
	float r = particle_distance(particle, other);
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
	return gaussian;
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
}
 
