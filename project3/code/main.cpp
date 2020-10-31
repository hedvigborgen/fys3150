#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
#include "time.h"
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

  double t;
  if (method == 1){
    clock_t start, finish;
    start = clock();

    Euler integrator(dt);
    solarSystem.writeToFile("../output/euler.xyz", "../output/euler_system.dat", "../output/angmom_euler.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      solarSystem.calculateAngMomentum();
      solarSystem.writeToFile("../output/euler.xyz", "../output/euler_system.dat", "../output/angmom_euler.dat", t);
      integrator.integrateOneStep(solarSystem);
    }

    finish = clock();
    double time = (double (finish - start)/CLOCKS_PER_SEC);
    cout << "Integration took " << time << " seconds to execute with Euler's method with n = " <<numTimesteps<< "." << endl;

  }

  else if (method == 2){
    clock_t start, finish;
    start = clock();

    VelocityVerlet integrator(dt);
    solarSystem.writeToFile("../output/verlet.xyz", "../output/verlet_system.dat", "../output/angmom_verlet.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      solarSystem.calculateAngMomentum();
      integrator.integrateOneStep(solarSystem);
      solarSystem.writeToFile("../output/verlet.xyz", "../output/verlet_system.dat", "../output/angmom_verlet.dat", t);

    }

    finish = clock();
    double time = (double (finish - start)/CLOCKS_PER_SEC);
    cout << "Integration took " << time << " seconds to execute with Verlet's method with n = " <<numTimesteps<< "." << endl;
  }

  //cout << "We just created a solar system that has " << solarSystem.bodies().size() << " objects." << endl;
  return 0;
}
