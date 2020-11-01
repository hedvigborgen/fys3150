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
    void calculateForcesAndEnergy(double beta);
    void calculateAngMomentum();
    // With *const we cannot set functin() = something.
    // The value can only be changed within a method
    int numberOfBodies() const;
    double totalEnergy() const;
    double potentialEnergy() const;
    double kineticEnergy() const;
    void writeToFile(std::string filename1, std::string filename2, std::string filename3, double t);
    void writeToFile_test(std::string filename, double t);
    std::vector<vec3> &angMomentum();
    std::vector<CelestialBody> &bodies();
    std::ofstream m_file_test;

private:
    // We cannot set m_value = something in main, since the value is private.
    // We can only reach the value through methods of the class.
    std::vector<vec3> m_angMomentum;
    std::vector<CelestialBody> m_bodies;
    std::ofstream m_file1, m_file2, m_file3;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double m_G;
};

#endif // SOLARSYSTEM_H

//Constructer == special method that is automatically called when an object of a class is created.
