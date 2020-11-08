#include "solarsystem.hpp"

// Constructor
SolarSystem::SolarSystem(){
  m_G = 4*M_PI*M_PI; // [AU^3/(M_sun yr^2)], Gravitational constant
}


// Creates celestial bodies
CelestialBody& SolarSystem::createCelestialBody(string name, vec3 position, vec3 velocity, double mass){
  m_bodies.push_back(CelestialBody(name, position, velocity, mass));
  return m_bodies.back();
}


// Returns number of celestial bodies
int SolarSystem::numberOfBodies() const{
  return m_bodies.size();
}


// Calculates position of the center of mass
void SolarSystem::calculateCenterofMass(){
  double tot_mass;
  for (int i=0; i<numberOfBodies(); i++){
    CelestialBody &body = m_bodies[i];
    m_posCenterofMass += body.mass*body.position;
    m_velCenterofMass += body.mass*body.velocity;
    tot_mass += body.mass;
  }
  m_posCenterofMass /= tot_mass;
  m_velCenterofMass /= tot_mass;
}


// Fixes the center of mass for the solar system
void SolarSystem::fixCenterofMass(){
  calculateCenterofMass();
  for(int i=0; i<numberOfBodies(); i++){
    CelestialBody &body = m_bodies[i];
    body.position -= m_posCenterofMass;
    body.velocity -= m_velCenterofMass;
  }
}


// Reads input file containing information about the celestial bodies
void SolarSystem::readinfo_SolarSystem(string fname, int numberOfBodies){
  vec3 pos;
  vec3 vel;
  double x, y, z, vx, vy, vz, mass;
  char name[numberOfBodies];
  const char *filename = fname.c_str();   //Each line of file gives initial condition for a particle on the form: x y z vx vy vz
  //Open files
  FILE *fp_init = fopen(filename, "r"); //Open file to read, specified by "r".

  for(int i=0; i < numberOfBodies; i++) {

    fscanf(fp_init, "%s %lf %lf %lf %lf %lf %lf %lf", name, &x, &y, &z, &vx, &vy, &vz, &mass); // One %lf (lf=long float or double) for each floating point number on each line of the file.
    pos = vec3(x, y, z);
    vel = vec3(vx*365.25, vy*365.25, vz*365.25);
    createCelestialBody(name, pos, vel, mass/1.989e30);
  }
  fixCenterofMass();
  fclose(fp_init); //Closing file with initial conditions
}


// Resets forces on celestial bodies and kinetic and potential energy for the system
void SolarSystem::resetForces_Energy(){
  m_kineticEnergy = 0.0;
  m_potentialEnergy = 0.0;
  for(CelestialBody &body : m_bodies){
    body.force.zeros();
  }
}


// Calculates forces acting on each celestial body, and the solar system's energy
void SolarSystem::calculateForcesAndEnergy(double beta){
  resetForces_Energy();
  vec3 deltaRVector;
  double M_1, M_2, dr;

  for (int i=0; i<numberOfBodies(); i++){
    CelestialBody &body1 = m_bodies[i];
    M_1 = body1.mass;

    for(int j=0; j<numberOfBodies(); j++){
      if (i != j){
        CelestialBody &body2 = m_bodies[j];
        M_2 = body2.mass;
        deltaRVector = body1.position - body2.position;
        dr = deltaRVector.length();
        body1.force -= m_G*M_1*M_2*deltaRVector/pow(dr, beta);
        m_potentialEnergy -= m_G*M_1*M_2/dr;
      }
    }
    m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
  }
}

// Calculates angular momentum for each celestial body
void SolarSystem::calculateAngMomentum(){
  m_angMomentum.clear();
  vec3 deltaRVector, pVector, angMom;
  double M_1, M_2, dr;

  for (int i=0; i<numberOfBodies(); i++){
    CelestialBody &body1 = m_bodies[i];
    M_1 = body1.mass;

    for(int j=0; j<numberOfBodies(); j++){

      if (i != j){
        CelestialBody &body2 = m_bodies[j];
        M_2 = body2.mass;
        deltaRVector = body1.position - body2.position;
        pVector = M_1*(body1.velocity - body2.velocity);
        angMom += deltaRVector.cross(pVector);
      }
    }
    m_angMomentum.push_back(angMom);
  }
}


double SolarSystem::totalEnergy() const{
  return m_kineticEnergy + m_potentialEnergy;
}


double SolarSystem::potentialEnergy() const{
  return m_potentialEnergy;
}


double SolarSystem::kineticEnergy() const{
  return m_kineticEnergy;
}

vector<vec3> SolarSystem::angularMomentum() const{
  return m_angMomentum;
}

vector<CelestialBody> &SolarSystem::bodies(){
  return m_bodies;
}
