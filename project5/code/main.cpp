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
  double charge, step, alpha, alpha0, deltaAlpha, omega;
  string write;

  dimension = 3;
  numberofParticles = 2;
  step = 1.0;


// Read in output file, abort if there are too few command-line arguments
if ((numArg != 7) && (numArg != 9)){
  cout << "Incorrect number of arguments!" << endl; exit(1);
}

else if (numArg == 7){
  maxVariations = atol(arguments[1]);
  whichMethod = atoi(arguments[2]);
  omega = atof(arguments[3]);
  equilibrationTime = atol(arguments[4]);
  MCCs = atol(arguments[5]);
  charge = atol(arguments[6]);

  if (whichMethod == 0){
   alpha = 1.0;
  }

  else if (whichMethod == 1){
    alpha = 0.85;
  }

  else if ((whichMethod != 0) && (whichMethod != 1)){
    cout << "Unacceptable choice of method!" << endl; exit(1);
  }
}

else if (numArg == 9){
  // Defining input arguments
  maxVariations = atol(arguments[1]);
  alpha0 = atof(arguments[2]);
  deltaAlpha = atof(arguments[3]);
  equilibrationTime = atol(arguments[4]);
  MCCs = atol(arguments[5]);
  charge = atol(arguments[6]);
  whichMethod = atoi(arguments[7]);
  write = arguments[8];
  omega = 1;
}


// Initializing the system with input arguments
QuantumDot quantumDot(dimension, numberofParticles, charge, equilibrationTime, MCCs, step);

if (maxVariations != 1){
   // Performing the MC sampling
  quantumDot.MonteCarlo(whichMethod, write, maxVariations, alpha0, deltaAlpha, omega);

  if (write == "at the end"){
    // Writing results to file after all MCCs
    quantumDot.WriteToFile(whichMethod);
  }
}

else if (maxVariations == 1){
  quantumDot.MonteCarlo(whichMethod, write, alpha, omega);
}

return 0;
}
