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

  SolarSystem solarSystem;
  MainFunc mainFunc;

  // Declaring input parameters
  int numTimesteps, numberOfBodies, method, choice;
  double dt;
  string fname = "";

  if (numArg == 7){
    choice = atoi(arguments[1]); // Determines if the code should;
    // 1: Use the normal gravitational force;
    // 2: Run for various betas, or;
    // 3: Compute the precession of Mercury
    numTimesteps = atoi(arguments[2]); // Number of timesteps
    dt = atof(arguments[3]); // Size of timestep

    if (choice != 3){
      fname = fname.append(arguments[4]); // Filename of inputfile
      numberOfBodies = atoi(arguments[5]); // Number of celestial bodies
      method = atoi(arguments[6]); // Euler or velocity Verlet method}

      solarSystem.readinfo_SolarSystem(fname, numberOfBodies);
      mainFunc.initializeBeta(choice);

      // Initializes the solar system
      mainFunc.initializeSolarSystem(solarSystem);

    }
  }


  // Asks for input arguments if not all were given
  else {
    cout << "Enter 1 for normal gravitational force," << endl;
    cout << "enter 2 to run script for various betas, or" << endl;
    cout << "enter 3 to calculate the precession of Mercury:" << endl;
    cin >> choice;
    cout << "Enter number of time steps:" << endl;
    cin >> numTimesteps;
    cout << "Enter size of time step:" << endl;
    cin >> dt;

    if (choice != 3){
      cout << "Enter input filename:" << endl;
      cin >> fname;
      cout << "Enter number of celestial bodies:" << endl;
      cin >> numberOfBodies;
      cout << "Enter 1 for forward Euler method, or enter 2 for velocity Verlet method:" << endl;
      cin >> method;

      solarSystem.readinfo_SolarSystem(fname, numberOfBodies);
      mainFunc.initializeBeta(choice);

      // Initializes the solar system
      mainFunc.initializeSolarSystem(solarSystem);

    }
  }


  // Executes with normal gravitational force for beta = 3.0
  if (choice == 1){

    // Calculates force and energies for first timestep
    solarSystem.calculateForcesAndEnergy(mainFunc.beta_vec[0], choice);

    clock_t start, finish; // Times the integration
    start = clock();

    // Loops over all timesteps and writes results to files
    mainFunc.timeLoop_reg(method, numTimesteps, dt, mainFunc.beta_vec[0], choice);

    finish = clock(); // Finish timing the integration
    double time = (double (finish - start)/CLOCKS_PER_SEC);
    cout << "Integration took " << time << " seconds to execute with n = " <<numTimesteps<< "." << endl;
  }


  // Tests gravitational force for different betas
  else if (choice == 2){

    // Loops over varying betas
    for (int idx = 0; idx < mainFunc.beta_vec.size(); idx++){

      // Calculates force and energies for first timestep
      solarSystem.calculateForcesAndEnergy(mainFunc.beta_vec[idx], choice);

      // Loops over all timesteps and writes results to files
      mainFunc.timeLoop_diffBeta(method, numTimesteps, dt, mainFunc.beta_vec[idx]);
    }
  }


  // Computes the precession of Mercury
  else if (choice == 3){
    fname = "../input/precession_mercury.txt";
    numberOfBodies = 2;
    method = 2; // Uses the velocity Verlet method

    solarSystem.readinfo_SolarSystem(fname, numberOfBodies);
    mainFunc.initializeBeta(1);

    // Initializes the solar system
    mainFunc.initializeSolarSystem(solarSystem);

    solarSystem.calculateForcesAndEnergy(mainFunc.beta_vec[0], choice);
    mainFunc.timeLoop_reg(method, numTimesteps, dt, mainFunc.beta_vec[0], choice);
  }

  return 0;
}
