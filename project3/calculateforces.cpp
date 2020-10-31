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
        }

        m_kineticEnergy += 0.5*body1.mass*body1.velocity.lengthSquared();
      }

    }
}