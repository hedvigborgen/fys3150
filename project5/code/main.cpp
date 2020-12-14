#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>


#include "QuantumDot.hpp"

using namespace std;
using namespace arma;

int main(int numArg, char *arguments[]){
  int dimension, numberofParticles, whichMethod;
  long int maxVariations, equilibrationTime, MCCs;
  double charge, step, alpha, alpha0, deltaAlpha, beta, omega;
  vec omegaVector;
  string write;

  dimension = 3;
  numberofParticles = 2;


// Read in output file, abort if there are too few command-line arguments
if ((numArg != 8) && (numArg != 10)){
  cout << "Incorrect number of arguments!" << endl; exit(1);
}

else if (numArg == 8){
  maxVariations = atol(arguments[1]);
  whichMethod = atoi(arguments[2]);
  equilibrationTime = atol(arguments[3]);
  MCCs = atol(arguments[4]);
  charge = atol(arguments[5]);
  step = atoi(arguments[6]);
  beta = atof(arguments[7]);

  if (whichMethod == 0){
   alpha = 1.0;
  }

  else if (whichMethod == 1 || whichMethod == 2){
    alpha = 0.85;
  }

  else if ((whichMethod != 0) && (whichMethod != 1)){
    cout << "Unacceptable choice of method!" << endl; exit(1);
  }
}


else if (numArg == 10){
  // Defining input arguments
  maxVariations = atol(arguments[1]);
  alpha0 = atof(arguments[2]);
  deltaAlpha = atof(arguments[3]);
  equilibrationTime = atol(arguments[4]);
  MCCs = atol(arguments[5]);
  step = atof(arguments[6]);
  charge = atol(arguments[7]);
  whichMethod = atoi(arguments[8]);
  write = arguments[9];
  omega = 1;
  beta = 1;
}


// Initializing the system with input arguments
QuantumDot quantumDot(dimension, numberofParticles, charge, equilibrationTime, MCCs, step, beta);

if (maxVariations != 1){
   // Performing the MC sampling
  quantumDot.MonteCarlo(whichMethod, write, maxVariations, alpha0, deltaAlpha, omega);

  if (write == "at the end"){
    // Writing results to file after all MCCs
    quantumDot.WriteToFile(whichMethod);
  }
}

else if (maxVariations == 1){
  double omegaList[] = {0.01, 0.5, 1.0};

  for (int i = 0; i < 3; i++){
    quantumDot.MonteCarlo(whichMethod, write, alpha, omegaList[i]);
  }
}




return 0;
}
