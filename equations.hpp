#ifndef _EQUATIONS_HPP
#define _EQUATIONS_HPP

#include "geometry.hpp"
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

/*float distance(Particle particle, Particle& other);
float gaussian_smoothing(Particle particle, Particle& other);
vec3 gradient_gaussian_smoothing(Particle particle, Particle& other);
float gradient2_gaussian_smoothing(Particle particle, Particle& other);
vec3 force_pressure(Particle particle, Particle& other, float particleMass);
vec3 force_viscosity(Particle particle, Particle& other, float particleMasss);*/
float particleDistance(Particle& current, Particle& other);
vec3 pressureForcePartial(Particle& current, Particle& other, float particleMass, float h);
vec3 viscosityForcePartial(Particle& current, Particle& other, float timeStepSize, float viscosityConstant, float particleMass, float h);
float spikyKernel(Particle& current, Particle& other, float h);
float spikyKernelGradient(Particle& current, Particle& other, float h);
float viscosityKernel(Particle& current, Particle& other, float h);
float viscosityKernelLaplacian(Particle& current, Particle& other, float h);
float poly6Kernel(Particle& current, Particle& other, float h);
vec3 poly6KernelGradient(Particle& current, Particle& other, float h);
float poly6KernelLaplacian(Particle& current, Particle& other, float h);
vec3 inwardSurfaceNormalPartial(Particle& current, Particle& other, float particleMass, float h);
float colorFieldLaplacianPartial(Particle& current, Particle& other, float particleMass, float h);

#endif
