#include "euler.hpp"
#include "solarsystem.hpp"

Euler::Euler(double dt)
{
  m_dt = dt;
}


void Euler::integrateOneStep(SolarSystem &system, double beta)
{   // Equivalent to "for body in bodies" in Python.
    for (CelestialBody &body : system.bodies()) {
        body.position += body.velocity*m_dt;
        body.velocity += body.force / body.mass * m_dt;
    }
    system.calculateForcesAndEnergy(beta);
}
