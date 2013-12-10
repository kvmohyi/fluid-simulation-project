#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#include <vector>

#include <fstream>
#include <sstream>
#include <iostream>

#include <glm/glm.hpp>

using namespace std;
using namespace glm;

class Particle {
public:
	float density;
	float pressure;
	vec3 position;
	vec3 velocity;
	vec3 acceleration;
    Particle(vec3 _position);
};

class Triangle {
public:
	vec3 a, b, c;
	vec3 normal;
        Triangle(vec3 a, vec3 b, vec3 c, vec3 normal);
};

class RigidBody {
public:
        float length, height, depth;
	vector<Triangle> triangles;
        RigidBody(float l, float h, float d);
        RigidBody();
	bool collision(vec3 start, vec3 end);
        pair<float, vec3> collisionTimeNormal(vec3 start, vec3 end);
        float rayTriangle(vec3 start, vec3 end, Triangle triangle);
};

#endif
