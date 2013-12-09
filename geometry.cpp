#include <vector>

#include <glm/glm.hpp>

#include "geometry.hpp"

using namespace std;
using namespace glm;

bool RigidBody::collision(vec3 start, vec3 end) {
	return false;
}

Particle::Particle(vec3 p){
  position = p;
}
