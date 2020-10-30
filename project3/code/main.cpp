#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
using namespace std;

int main(int numArguments, char **arguments)
{
    int numTimesteps = 1000;
    if(numArguments >= 2) numTimesteps = atoi(arguments[1]);

    SolarSystem solarSystem;
    // To store the referance: CelestialBody &sun = solarSystem.createCelestialBody( vec3, vec3, mass);
    solarSystem.createCelestialBody( vec3(0,0,0), vec3(0,0,0), 1.0 ); //Sun
    solarSystem.createCelestialBody( vec3(1, 0, 0), vec3(0, 2*M_PI, 0), 3e-6 ); //Earth

    // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
    vector<CelestialBody> &bodies = solarSystem.bodies();

    for(int i = 0; i<bodies.size(); i++) {
        CelestialBody &body = bodies[i]; // Reference to this body
        cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    }

    double dt = 0.001;
    // Euler integrator(dt);
    VelocityVerlet integrator(dt);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
        // For the integrator to do as we want we need to calculate the forces in between the timesteps
        solarSystem.calculateForcesAndEnergy();
        integrator.integrateOneStep(solarSystem);
        // solarSystem.writeToFile("../output/euler.xyz");
        solarSystem.writeToFile("../output/verlet.xyz");
    }

    cout << "We just created a solar system that has " << solarSystem.bodies().size() << " objects." << endl;
    return 0;
}