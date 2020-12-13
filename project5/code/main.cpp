#include <iostream>
#include <fstream>
#include <iomanip>
// #include "lib.h"

#include "QuantumDot.hpp"

using namespace std;
using namespace arma;

int main(int numArg, char *arguments[]){
  int dimension, numberofParticles, whichMethod;
  long int maxVariations, equilibrationTime, MCCs;
  double charge, step, alpha0, deltaAlpha, omega;
  string write;

  dimension = 3;
  numberofParticles = 2;
  step = 1.0;
  // alpha0 = 0.6;
  // deltaAlpha = 0.2;
  omega = 1;


// Read in output file, abort if there are too few command-line arguments
if (numArg != 9){
  cout << "Incorrect number of arguments!" << endl; exit(1);
}

else {
  // Defining input arguments
  maxVariations = atol(arguments[1]);
  alpha0 = atof(arguments[2]);
  deltaAlpha = atof(arguments[3]);
  equilibrationTime = atol(arguments[4]);
  MCCs = atol(arguments[5]);
  charge = atol(arguments[6]);
  whichMethod = atoi(arguments[7]);
  write = arguments[8];
}

// Initializing the system, read in data
QuantumDot quantumDot(dimension, numberofParticles, charge, maxVariations,
  equilibrationTime, MCCs, step, alpha0, deltaAlpha, omega);

quantumDot.Initialize();

 // Do the mc sampling
quantumDot.MonteCarlo(whichMethod, write);

// Print out results
quantumDot.WriteToFile(whichMethod);


return 0;

}
