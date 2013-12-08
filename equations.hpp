#ifndef _EQUATIONS_HPP
#define _EQUATIONS_HPP

#include "simulation.hpp"
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

float distance(Particle particle, Particle other);
float gaussian_smoothing(Particle particle, Particle other);
vec3 gradient_gaussian_smoothing(Particle particle, Particle other);
float gradient2_gaussian_smoothing(Particle particle, Particle other);
vec3 force_pressure(Particle particle, Particle other);
vec3 force_viscosity(Particle particle, Particle other);

#endif
