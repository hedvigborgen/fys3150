#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include "celestialbody.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

class SolarSystem
{
public:
    SolarSystem(); // Constructer. Equivalent to __init__() in Python.
    CelestialBody &createCelestialBody(vec3 position, vec3 velocity, double mass);
    void calculateForcesAndEnergy();
    // With *const we cannot set functin() = something.
    // The value can only be changed within a method
    int numberOfBodies() const;
    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename);
    vec3 angularMomentum() const;
    // Vector, named bodies, with different celestial bodies as elements.
    std::vector<CelestialBody> &bodies();

private:
    // We cannot set m_value = something in main, since the value is private.
    // We can only reach the value through methods of the class.
    std::vector<CelestialBody> m_bodies;
    vec3 m_angularMomentum;
    std::ofstream m_file;
    double m_kineticEnergy;
    double m_potentialEnergy;
};

#endif // SOLARSYSTEM_H

//Constructer == special method that is automatically called when an object of a class is created.
