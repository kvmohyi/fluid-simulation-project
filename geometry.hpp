#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <vector>

#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Particle {
public:
	float mass;
	float massDensity;
	float pressure;
	vec3 position;
	vec3 velocity;
	vec3 acceleration;
};

class Triangle {
public:
	vec3 a, b, c, d;
	vec3 normal;
};

class RigidBody {
public:
	vector<Triangle> triangles;
	vec3 normal;

	bool collision(vec3 start, vec3 end);
};

#endif