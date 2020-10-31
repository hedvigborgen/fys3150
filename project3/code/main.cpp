#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
using namespace std;

int main(int numArg, char **arguments)
{
  int numTimesteps;
  double dt;
  int method;
  int arg;

  if (numArg == 4) {
    numTimesteps = atoi(arguments[1]);
    dt = atof(arguments[2]);
    method = atoi(arguments[3]);
  } else {
    cout << "Enter number of time steps:" << endl;
    cin >> numTimesteps;
    cout << "Enter value for time step:" << endl;
    cin >> dt;
    cout << "Enter 1 for forward Euler method, or enter 2 for velocity Verlet method:" << endl;
    cin >> arg;
  }

  SolarSystem solarSystem;
  // To store the referance: CelestialBody &sun = solarSystem.createCelestialBody( vec3, vec3, mass);
  solarSystem.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1.0); //Sun
  solarSystem.createCelestialBody(vec3(1,0,0), vec3(0,2*M_PI,0), 3e-6); //Earth
  solarSystem.calculateForcesAndEnergy();

  // To get a list (a reference, not copy) of all the bodies in the solar system, we use the .bodies() function
  vector<CelestialBody> &bodies = solarSystem.bodies();

    // To check that the positions and velocities of the bodies are initialized
    //for(int i = 0; i<bodies.size(); i++) {
        //CelestialBody &body = bodies[i]; // Reference to this body
        //cout << "The position of this object is " << body.position << " with velocity " << body.velocity << endl;
    //}


  double t;
  if (method == 1){
    // clock_t start, finish;
    // start = clock();

    Euler integrator(dt);
    solarSystem.writeToFile("../output/euler.xyz", "../output/euler_system.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      solarSystem.writeToFile("../output/euler.xyz", "../output/euler_system.dat", t);
      integrator.integrateOneStep(solarSystem);
    }
    // finish = clock();
    // double time = (double (finish - start)/CLOCKS_PER_SEC);

  } else if (method == 2){
    VelocityVerlet integrator(dt);
    solarSystem.writeToFile("../output/verlet.xyz", "../output/verlet_system.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      integrator.integrateOneStep(solarSystem);
      solarSystem.writeToFile("../output/verlet.xyz", "../output/verlet_system.dat", t);
    }
  }

  //cout << "We just created a solar system that has " << solarSystem.bodies().size() << " objects." << endl;
  return 0;
}
