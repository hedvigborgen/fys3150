#include "solarsystem.hpp"
#include <iostream>
using namespace std;

SolarSystem::SolarSystem()    // Constructer. Equivalent to __init__() in Python.
{
  m_kineticEnergy = 0;
  m_potentialEnergy = 0;
}

CelestialBody& SolarSystem::createCelestialBody(vec3 position, vec3 velocity, double mass) {
    m_bodies.push_back( CelestialBody(position, velocity, mass) ); // push_back() = append() in Python
    return m_bodies.back(); // Return reference to the newest added celstial body
}

void SolarSystem::calculateForcesAndEnergy()
{
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
    m_angularMomentum.zeros();
    double G = 4*M_PI*M_PI; // [AU^3/(M_sun yr^2)], Gravitational constant

    for(CelestialBody &body : m_bodies) {
        // Reset forces on all bodies
        body.force.zeros();
    }
    vec3 deltaRVector;
    double M_1, M_2, dr;

    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body1 = m_bodies[i];
        M_1 = body1.mass;
        for(int j=0; j<numberOfBodies(); j++) {
          if (i != j){
            CelestialBody &body2 = m_bodies[j];
            M_2 = body2.mass;
            deltaRVector = body1.position - body2.position;
            dr = deltaRVector.length();
            body1.force -= G*M_1*M_2*deltaRVector/(dr*dr*dr);
            m_potentialEnergy -= G*M_1*M_2/dr;
            // body1.m_potentialEnergy -= G*M_1*M_2/dr;
            // body2.m_potentialEnergy += body1.m_potentialEnergy;
        }
        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
      }
    }
}

int SolarSystem::numberOfBodies() const
{
    return m_bodies.size();
}

double SolarSystem::totalEnergy() const
{
    return m_kineticEnergy + m_potentialEnergy;
}

double SolarSystem::potentialEnergy() const
{
    return m_potentialEnergy;
}

double SolarSystem::kineticEnergy() const
{
    return m_kineticEnergy;
}

void SolarSystem::writeToFile(string filename1, string filename2, double t)
{
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

    for(CelestialBody &body : m_bodies) {
        m_file1 << t << " "
          << body.position.x() << " "
          << body.position.y() << " "
          << body.position.z() << "\n";
    }
    m_file2 << t << " " << m_potentialEnergy << " " << m_kineticEnergy << " " << totalEnergy() << "\n";
}

vec3 SolarSystem::angularMomentum() const
{
    return m_angularMomentum;
}

std::vector<CelestialBody> &SolarSystem::bodies()
{
    return m_bodies;
}