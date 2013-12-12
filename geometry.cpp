#include <vector>
#include <utility> 

#include <glm/glm.hpp>

#include "geometry.hpp"

using namespace std;
using namespace glm;

Triangle::Triangle(vec3 va, vec3 vb, vec3 vc, vec3 n){
  a = va;
  b = vb;
  c = vc;
  normal = n;
}
RigidBody::RigidBody(){
}

RigidBody::RigidBody(float l, float h, float d){
  length = l;
  height = h;
  depth = d;
  float halfLength = length / 2.0;
  float halfHeight = height / 2.0;
  float halfDepth = depth / 2.0;
  vec3 vertex1 = vec3(halfLength, halfHeight, halfDepth);
  vec3 vertex2 = vec3(-1 * halfLength, halfHeight, halfDepth);
  vec3 vertex3 = vec3(-1 * halfLength, halfHeight, -1 * halfDepth);
  vec3 vertex4 = vec3(halfLength, halfHeight, -1 * halfDepth);
  vec3 vertex5 = vec3(halfLength, -1 * halfHeight, halfDepth);
  vec3 vertex6 = vec3(-1 * halfLength, -1 * halfHeight, halfDepth);
  vec3 vertex7 = vec3(-1 * halfLength, -1 * halfHeight, -1 * halfDepth);
  vec3 vertex8 = vec3(halfLength, -1 * halfHeight, -1 * halfDepth);
  triangles.push_back(Triangle(vertex6, vertex1, vertex5, vec3(0.0, 0.0, -1.0)));
  triangles.push_back(Triangle(vertex6, vertex2, vertex1, vec3(0.0, 0.0, -1.0)));
  triangles.push_back(Triangle(vertex5, vertex4, vertex8, vec3(-1.0, 0.0, 0.0)));
  triangles.push_back(Triangle(vertex5, vertex1, vertex4, vec3(-1.0, 0.0, 0.0)));
  triangles.push_back(Triangle(vertex2, vertex4, vertex1, vec3(0.0, -1.0, 0.0)));
  triangles.push_back(Triangle(vertex2, vertex3, vertex4, vec3(0.0, -1.0, 0.0)));
  triangles.push_back(Triangle(vertex6, vertex7, vertex3, vec3(1.0, 0.0, 0.0)));
  triangles.push_back(Triangle(vertex6, vertex3, vertex2, vec3(1.0, 0.0, 0.0)));
  triangles.push_back(Triangle(vertex7, vertex8, vertex4, vec3(0.0, 0.0, 1.0)));
  triangles.push_back(Triangle(vertex7, vertex4, vertex3, vec3(0.0, 0.0, 1.0)));
  triangles.push_back(Triangle(vertex6, vertex5, vertex8, vec3(0.0, 1.0, 0.0)));
  triangles.push_back(Triangle(vertex6, vertex8, vertex7, vec3(0.0, 1.0, 0.0)));
}

void RigidBody::handleCollisions(Particle& particle) {
  float dampening = 0.8f;
  if (particle.position.x < -length / 2.0) {
    particle.position.x = -length / 2.0 + (-length / 2.0 - particle.position.x);
    particle.velocity.x = -1.0 * dampening * particle.velocity.x;
  }
  else if (particle.position.x > length / 2.0) {
    particle.position.x = length / 2.0 - (-length / 2.0 + particle.position.x);
    particle.velocity.x = -1.0 * dampening * particle.velocity.x;
  }

  if (particle.position.y < -length / 2.0) {
    particle.position.y = -length / 2.0 + (-length / 2.0 - particle.position.y);
    particle.velocity.y = -1.0 * dampening * particle.velocity.y;
  }
  else if (particle.position.y > length / 2.0) {
    particle.position.y = length / 2.0 - (-length / 2.0 + particle.position.y);
    particle.velocity.y = -1.0 * dampening * particle.velocity.y;
  }

  if (particle.position.z < -length / 2.0) {
    particle.position.z = -length / 2.0 + (-length / 2.0 - particle.position.z);
    particle.velocity.z = -1.0 * dampening * particle.velocity.z;
  }
  else if (particle.position.z > length / 2.0) {
    particle.position.z = length / 2.0 - (-length / 2.0 + particle.position.z);
    particle.velocity.z = -1.0 * dampening * particle.velocity.z;
  }
}

bool RigidBody::collision(vec3 start, vec3 end) {
  if(end.x > length / 2.0 || end.x < -1 * length / 2.0){
    return true;
  }
  else if(end.y > height / 2.0 || end.y < -1 * height / 2.0){
    return true;
  }
  else if(end.z > depth / 2.0 || end.z < -1 * depth / 2.0){
    return true;
  }
  return false;
}

pair<float, vec3> RigidBody::collisionTimeNormal(vec3 start, vec3 end) {
  for (int s = 0; s < triangles.size(); s++){
    float time = rayTriangle(start, end, triangles[s]);
    if(time >= 0.0){
      vec3 normal = triangles[s].normal;
      pair <float, vec3> toReturn = make_pair(time, normal);
      return toReturn;
    }
  }
  return make_pair(0.0, vec3(1.0, 0.0, 0.0));
}

float RigidBody::rayTriangle(vec3 start, vec3 end, Triangle triangle){
  vec3 rayStep = end - start;
  float a = triangle.a.x - triangle.b.x;
  float b = triangle.a.y - triangle.b.y;
  float c = triangle.a.z - triangle.b.z;
  float d = triangle.a.x - triangle.c.x;
  float e = triangle.a.y - triangle.c.y;
  float f = triangle.a.z - triangle.c.z;
  float g = rayStep.x;
  float h = rayStep.y;
  float i = rayStep.z;
  float j = triangle.a.x - start.x;
  float k = triangle.a.y - start.y;
  float l = triangle.a.z - start.z;
  float m = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);
  float time = (f * (a * k - j*b) + e * (j * c - a * l) + d * (b * l - k * c)) / m;
  time = time * -1;
  if (time < 0.0 || time > 1.0) {
    return -1.0;
  }
  return time;
}								
Particle::Particle(vec3 p){
  position = p;
}

Particle::Particle(vec3 _position, vec3 _velocity) {
  position = _position;
  velocity = _velocity;
}

Particle::Particle(vec3 _position, vec3 _velocity, float _density, float _pressure) {
  position = _position;
  velocity = _velocity;
  density = _density;
  pressure = _pressure;
}

void Particle::print() {
  cout << "Particle" << endl;
  cout << "\tPosition: " << position.x << ", " << position.y << ", " << position.z << endl;
  cout << "\tvelocity: " << velocity.x << ", " << velocity.y << ", " << velocity.z << endl;
  cout << "\tPressure: " << pressure << endl;
  cout << "\tDensity: " << density << endl;
}
