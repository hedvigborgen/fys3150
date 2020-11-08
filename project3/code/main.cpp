#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
#include "mainfunc.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "time.h"
using namespace std;


int main(int numArg, char *arguments[]){
  int numTimesteps, numberOfBodies, method, choice;
  double dt;
  string fname = "";

  if (numArg == 7){
    numTimesteps = atoi(arguments[1]);
    dt = atof(arguments[2]);
    fname = fname.append(arguments[3]);
    numberOfBodies = atoi(arguments[4]);
    method = atoi(arguments[5]);
    choice = atoi(arguments[6]);
  }

  else {
    cout << "Enter number of time steps:" << endl;
    cin >> numTimesteps;
    cout << "Enter size of time step:" << endl;
    cin >> dt;
    cout << "Enter input filename:" << endl;
    cin >> fname;
    cout << "Enter number of celestial bodies:" << endl;
    cin >> numberOfBodies;
    cout << "Enter 1 for forward Euler method, or enter 2 for velocity Verlet method:" << endl;
    cin >> method;
    cout << "Enter 1 for normal gravitational force or 2 for test of different forces:" << endl;
    cin >> choice;
  }


SolarSystem solarSystem;
MainFunc mainFunc;

solarSystem.readinfo_SolarSystem(fname, numberOfBodies);
mainFunc.initializeBeta(choice);


// Executing with normal gravitational force for beta = 3.0
if (choice == 1){
  // Initializing
  solarSystem.calculateForcesAndEnergy(mainFunc.beta_vec[0]);

// Timing the integration
  clock_t start, finish;
  start = clock();

// Integration
  mainFunc.integration1(method, numTimesteps, dt, mainFunc.beta_vec[0]);

  finish = clock();
  double time = (double (finish - start)/CLOCKS_PER_SEC);
  cout << "Integration took " << time << " seconds to execute with Verlet's method with n = " <<numTimesteps<< "." << endl;
}


// Testing gravitational force for different betas
else if (choice == 2){
  for (int idx = 0; idx < mainFunc.beta_vec.size(); idx++){
    solarSystem.calculateForcesAndEnergy(mainFunc.beta_vec[idx]);
    mainFunc.integration2(method, numTimesteps, dt, mainFunc.beta_vec[idx]);
  }
}

return 0;
}
