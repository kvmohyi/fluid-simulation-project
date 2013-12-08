
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include <glm/glm.hpp>

#include 'simulations.cpp'

#define PI 3.14159265  

using namespace std;
using namespace glm;


float sqr(float x);
float distance(Particle particle, Particle other);
float gaussian_smoothing(Particle particle, Particle other);
vec3 gradient_gaussian_smoothing(Particle particle, Particle other);
vec3 gradient_gaussian_smoothing(Particle particle, Particle other);
float density(Particle particle);
vec3 pressure(Particle particle);
vec3 force_pressure(Particle particle, Particle other);
vec3 force_viscosity(Particle particle, Particle other);

float sqr(float x){
  return x * x;
}

float distance(Particle particle, Particle other){
        vec3 pos1 = particle.position;
	vec3 pos2 = other.position;
	return sqrt(sqr(pos1.x - pos2.x) + sqr(pos1.y - pos2.y) + sqr(pos1.z - pos2.z));
}		    
float gaussian_smoothing(Particle particle, Particle other){
	float r = distance(particle, other);
	if r > radius {
		return 0;
	}
	float gaussian = exp(-4 * sqr(r/(float)(radius)));
	return gaussian;
}

vec3 gradient_gaussian_smoothing(Particle particle, Particle other){
	float r = distance(particle, other);
	if r > radius {
		return vec3(0.0, 0.0, 0.0);
	}
	float gaussian = gaussian_smoothing(particle, other);  
	float partialX = -8 * abs(particle.position.x - other.position.x);
	float partialY = -8 * abs(particle.position.y - other.position.y);
	float partialZ = -8 * abs(particle.position.z - other.position.z);
	vec3 gradient = new vec3(partialX * gaussian, partialY * gaussian, partialZ * gaussian);
	return gradient;

}

vec3 gradient2_gaussian_smoothing(Particle particle, Particle other){	
	float r = distance(particle, other);
	if r > radius {
		return vec3(0.0, 0.0, 0.0);
	}
	float gaussian = gaussian_smoothing(particle, other);
	float partialX = -8 * gaussian + 64 * sqr(abs(particle.position.x - other.position.x)) * gaussian;
	float partialY = -8 * gaussian + 64 * sqr(abs(particle.position.y - other.position.y)) * gaussian;
	float partialZ = -8 * gaussian + 64 * sqr(abs(particle.position.z - other.position.z)) * gaussian;
	vec3 gradient2 = new vec3(partialX, partialY, partialZ);
	
}
	
vec3 force_pressure(Particle particle, Particle other) {
	float avgPressure = (particle.pressure + other.pressure) / 2.0;
	avgPressure = avgPressure * other.mass / other.mass_density;
	vec3 force = gradient_gaussian_smoothing(particle, other);
	force = force * avgPressure;
	return force;
}

vec3 force_viscosity(Particle particle, Particle other) {
	float diffVelocity = other.velocity - particle.velocity;
	diffVelocity = diffVelocity * other.mass / other.mass_density;
	vec3 force = gradient2_gaussian_smoothing(particle, other);
	force = force * diffVelocity;
	return force;
}
 
