#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H

#include "vec3.hpp"
using namespace std;

class CelestialBody {
public:
  string nameOfBody;
  vec3 position;
  vec3 velocity;
  vec3 force;
  double mass;

  CelestialBody(string name, vec3 position, vec3 velocity, double mass);
  void resetForce();
};

#endif // CELESTIALBODY_H
