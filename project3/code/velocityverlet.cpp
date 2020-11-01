#include "velocityverlet.hpp"
#include "solarsystem.hpp"

VelocityVerlet::VelocityVerlet(double dt)
{
  m_dt = dt;
}

void VelocityVerlet::integrateOneStep(SolarSystem &system, double beta)
{
    for(CelestialBody &body : system.bodies()) {
        // update position to x_i+1
        body.position += m_dt*body.velocity+m_dt*m_dt/2*(body.force / body.mass);
        // v_i+1 = v_i + a_i*h/2
        body.velocity += body.force / body.mass * m_dt/2;
    }
    system.calculateForcesAndEnergy(beta);
    for (CelestialBody &body : system.bodies()) {
        // v_i+1 += a_i+1*h/2
        body.velocity += body.force / body.mass * m_dt/2;
    }
}
