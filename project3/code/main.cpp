#include <iostream>
#include <cmath>
#include <cstdlib>
#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
#include "time.h"
using namespace std;

int main(int numArg, char **arguments){
  int numTimesteps, method, choice;
  double dt;


  if (numArg == 5){
    numTimesteps = atoi(arguments[1]);
    dt = atof(arguments[2]);
    method = atoi(arguments[3]);
    choice = atoi(arguments[4]);
  }

  else {
    cout << "Enter number of time steps:" << endl;
    cin >> numTimesteps;
    cout << "Enter value for time step:" << endl;
    cin >> dt;
    cout << "Enter 1 for forward Euler method, or enter 2 for velocity Verlet method:" << endl;
    cin >> method;
    cout << "Enter 1 for normal gravitational force or 2 for test of different forces:" << endl;
    cin >> choice;
  }

SolarSystem solarSystem;
vector<CelestialBody> &bodies = solarSystem.bodies();

double t;

if (choice == 1){
  double beta = 3.0;

  // To store the referance: CelestialBody &sun = solarSystem.createCelestialBody( vec3, vec3, mass);
  solarSystem.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1.0); //Sun
  solarSystem.createCelestialBody(vec3(1,0,0), vec3(0,2*M_PI,0), 3e-6); //Earth
  solarSystem.calculateForcesAndEnergy(beta);

  if (method == 1){
    clock_t start, finish;
    start = clock();

    Euler integrator(dt);
    solarSystem.writeToFile("../output/euler_positions.xyz", "../output/euler_energies.dat", "../output/euler_angmom.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      solarSystem.calculateAngMomentum();
      solarSystem.writeToFile("../output/euler_positions.xyz", "../output/euler_energies.dat", "../output/euler_angmom.dat", t);
      integrator.integrateOneStep(solarSystem, beta);
    }

    finish = clock();
    double time = (double (finish - start)/CLOCKS_PER_SEC);
    cout << "Integration took " << time << " seconds to execute with Euler's method with n = " <<numTimesteps<< "." << endl;

  }

  else if (method == 2){
    clock_t start, finish;
    start = clock();

    VelocityVerlet integrator(dt);
    solarSystem.writeToFile("../output/verlet_positions.xyz", "../output/verlet_energies.dat", "../output/verlet_angmom.dat", 0);
    for(int timestep=0; timestep<numTimesteps; timestep++) {
      t = timestep*dt;
      solarSystem.calculateAngMomentum();
      integrator.integrateOneStep(solarSystem, beta);
      solarSystem.writeToFile("../output/verlet_positions.xyz", "../output/verlet_energies.dat", "../output/verlet_angmom.dat", t);

    }

    finish = clock();
    double time = (double (finish - start)/CLOCKS_PER_SEC);
    cout << "Integration took " << time << " seconds to execute with Verlet's method with n = " <<numTimesteps<< "." << endl;
  }
}




// Testing force for different betas
else if (choice == 2){
  int idx = 0;
  double beta[6];
  for (int b=20; b<31; b+=2){
    beta[idx] = b/10.0;
    idx += 1;
  }

  for (int idx = 0; idx < 6; idx++){
    // To store the referance: CelestialBody &sun = solarSystem.createCelestialBody( vec3, vec3, mass);
    solarSystem.createCelestialBody(vec3(0,0,0), vec3(0,0,0), 1.0); //Sun
    solarSystem.createCelestialBody(vec3(1,0,0), vec3(0,2*M_PI,0), 3e-6); //Earth
    solarSystem.calculateForcesAndEnergy(beta[idx]);

    if (method == 1){
      Euler integrator(dt);
      solarSystem.m_file_test << beta[idx] << endl;
      solarSystem.writeToFile_test("../output/euler_test_positions_beta.xyz", 0);
      for(int timestep=0; timestep<numTimesteps; timestep++) {
        t = timestep*dt;
        solarSystem.calculateAngMomentum();
        solarSystem.writeToFile_test("../output/euler_test_positions_beta.xyz", t);
        integrator.integrateOneStep(solarSystem, beta[idx]);
      }
    }

    else if (method == 2){
      solarSystem.m_file_test << beta[idx] << endl;
      VelocityVerlet integrator(dt);
      solarSystem.writeToFile_test("../output/verlet_test_positions_beta.xyz", 0);
      for(int timestep=0; timestep<numTimesteps; timestep++) {
        t = timestep*dt;
        solarSystem.calculateAngMomentum();
        integrator.integrateOneStep(solarSystem, beta[idx]);
        solarSystem.writeToFile_test("../output/verlet_test_positions_beta.xyz", t);
      }
    }
  }
}

return 0;
}
//linje 106, 110, 118, 123 må fikses
