#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
using namespace std;


class SolarSystem {
public:
  SolarSystem();
  CelestialBody &createCelestialBody(string name, vec3 position, vec3 velocity, double mass);
  void readinfo_SolarSystem(string fname, int numberOfBodies);
  void calculateForcesAndEnergy(double beta, int choice);
  void calculateAngMomentum();
  int numberOfBodies() const;
  double totalEnergy() const;
  double potentialEnergy() const;
  double kineticEnergy() const;
  vector<vec3> angularMomentum() const;
  vector<CelestialBody> &bodies();


private:
  vector<CelestialBody> m_bodies;
  vec3 m_posCenterofMass;
  vec3 m_velCenterofMass;
  vec3 m_angMom;
  vector<vec3> m_angMomentum;
  double m_kineticEnergy;
  double m_potentialEnergy;
  double m_G;
  void calculateCenterofMass();
  void fixCenterofMass();
  void resetForces_Energy();
};

#endif // SOLARSYSTEM_H
