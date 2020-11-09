#include "celestialbody.hpp"


CelestialBody::CelestialBody(string name, vec3 pos, vec3 vel, double mass_){
  nameOfBody = name;
  position = pos;
  velocity = vel;
  mass = mass_;
}

void CelestialBody::resetForce(){
    force.zeros();
}
