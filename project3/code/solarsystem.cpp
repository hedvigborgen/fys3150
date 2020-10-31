#include "solarsystem.hpp"
#include <iostream>
using namespace std;

SolarSystem::SolarSystem(){
  m_G = 4*M_PI*M_PI; // [AU^3/(M_sun yr^2)], Gravitational constant
  m_kineticEnergy = 0;
  m_potentialEnergy = 0;
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass){
    m_bodies.push_back( CelestialBody(position, velocity, mass) ); // push_back() = append() in Python
    return m_bodies.back(); // Return reference to the newest added celestial body
}

void SolarSystem::calculateForcesAndEnergy(){
    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }
    vec3 deltaRVector;
    double M_1, M_2, dr;
    int idx = 0;
    double beta = new double[n];
    for (int i=20; i<=30; i++){
      beta[idx] = i/10
      idx += 1
    }

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        M_1 = body1.mass;

        for(int j=0; j<numberOfBodies(); j++) {

          if (i != j){
            CelestialBody &body2 = m_bodies[j];
            M_2 = body2.mass;
            deltaRVector = body1.position - body2.position;
            dr = deltaRVector.length();
            body1.force -= m_G*M_1*M_2*deltaRVector/pow(dr,beta);
            m_potentialEnergy -= m_G*M_1*M_2/dr;
          }
        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
      }
    }
}


void SolarSystem::calculateAngMomentum(){

  double M_sun, M_planet, dr;
  CelestialBody &body_S = m_bodies[0];
  M_sun = body_S.mass;
  m_angMomentum.clear();
  for(int i=1; i<numberOfBodies(); i++) {
    CelestialBody &body_p = m_bodies[i];
    M_planet = body_p.mass;
    vec3 r = body_S.position - body_p.position;
    vec3 p = M_planet*body_p.velocity;
    m_angMomentum.push_back(r.cross(p));
  }
}


int SolarSystem::numberOfBodies() const{
    return m_bodies.size();
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

void SolarSystem::writeToFile(string filename1, string filename2, string filename3, double t){
    if (!m_file1.good()) {
        m_file1.open(filename1.c_str(), ofstream::out);
        if(!m_file1.good()) {
            cout << "Error opening file " << filename1 << ". Aborting!" << endl;
            terminate();
        }
    }
    if (!m_file2.good()) {
        m_file2.open(filename2.c_str(), ofstream::out);
        if(!m_file2.good()) {
            cout << "Error opening file " << filename2 << ". Aborting!" << endl;
            terminate();
        }
    }
    if (!m_file3.good()) {
        m_file3.open(filename3.c_str(), ofstream::out);
        if(!m_file3.good()) {
            cout << "Error opening file " << filename3 << ". Aborting!" << endl;
            terminate();
        }
    }

    for(CelestialBody &body : m_bodies) {
        m_file1 << t << " "
          << body.position.x() << " "
          << body.position.y() << " "
          << body.position.z() << "\n";
    }

    m_file2 << t << " " << m_potentialEnergy << " " << m_kineticEnergy << " " << totalEnergy() << "\n";

    for(vec3 &angMom : m_angMomentum){
      m_file3 << t << " "
        << angMom(0) << " "
        << angMom(1) << " "
        << angMom(2) << "\n";
    }
}
std::vector<CelestialBody> &SolarSystem::bodies(){
    return m_bodies;
}
