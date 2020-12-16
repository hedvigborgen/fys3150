#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <armadillo>
#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "QuantumDot.hpp"

using namespace std;
using namespace arma;


int main(int numArg, char *arguments[]){
  int dimension, numberofParticles, whichMethod;
  long int MCCs, equilibrationTime, maxVariations;
  double step, charge, alpha, alpha0, deltaAlpha, beta, omega, omega0, deltaOmega;
  string task, write;

  // Defining set parameters for all tasks
  dimension = 3;
  numberofParticles = 2;
  task = arguments[1];

  // For studying the stability of the algorithm as function of Monte Carlo cycles
  if (task == "MCCs"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    maxVariations = atol(arguments[4]);
    alpha0 = atof(arguments[5]);
    deltaAlpha = atof(arguments[6]);

    // Predecided input arguments
    equilibrationTime = 0;
    write = "each time";
    step = 1.0;
    beta = 1.0;
    omega = 1.0;
    charge = 1.0;
  }

  // For deciding the optimal value of the step length as function of alpha
  else if (task == "StepLength"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    maxVariations = atol(arguments[4]);
    alpha0 = atof(arguments[5]);
    deltaAlpha = atof(arguments[6]);
    step = atof(arguments[7]);

    // Predecided input arguments
    equilibrationTime = 100000;
    write = "at the end";
    beta = 1.0;
    omega = 1.0;
    charge = 1.0;
  }

  // For finding the expectation value of the energy and the variance
  // for certain set values of the parameters alpha, beta and omega
  else if (task == "Parameters"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    alpha = atof(arguments[4]);
    beta = atof(arguments[5]);
    omega = atof(arguments[6]);

    // Predecided input arguments
    maxVariations = 1;
    equilibrationTime = 100000;
    charge = 1.0;

    if ((whichMethod != 0) && (whichMethod != 1) && (whichMethod != 2)){
      cout << "Unacceptable choice of method!" << endl; exit(1);
    }
  }

  // For finding the expectation value of the energy and the variance
  // as function of the parameters alpha and beta with one certain value of omega
  else if (task == "Loop"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega = atof(arguments[3]);

    // Predecided input arguments
    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 2;
    alpha0 = 0.60;
    deltaAlpha = 0.01;
    charge = 1.0;
  }

  // For testing the compliance of the energies calculated without Coulomb interaction
  // with the Virial theorem
  else if (task == "VirialwithoutInteraction"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega0 = atof(arguments[3]);
    deltaOmega = atof(arguments[4]);

    // Predecided input arguments
    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 0;
    alpha = 0.995;
    beta = 0.280;
    charge = 1.0;
  }

  // For testing the compliance of the energies calculated with Coulomb interaction
  // with the Virial theorem
  else if (task == "VirialwithInteraction"){

    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega0 = atof(arguments[3]);
    deltaOmega = atof(arguments[4]);

    // Predecided input arguments
    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 2;
    alpha = 0.995;
    beta = 0.280;
    charge = 1.0;
  }


  // Quitting the program if unvalid entries for "task" is provided
  else if ((task != "MCCs") && (task != "StepLength") && (task != "Parameters")
  && (task != "Loop") && (task != "VirialwithoutInteraction") && (task != "VirialwithInteraction")){
    cout << "Unacceptable choice of task!" << endl; exit(1);
  }


  // Initializing the quantum system with input arguments
  QuantumDot quantumDot(dimension, numberofParticles, charge, equilibrationTime, MCCs);


  // Performing the MC sampling for different values of the parameter alpha
  // and writing output files
  if (maxVariations != 1){
    quantumDot.MonteCarlo(whichMethod, write, maxVariations, alpha0, deltaAlpha, beta, omega, step);
  }


  // Performing the MC sampling for certain values of the parameters alpha,
  // beta and omega
  else if (task == "Parameters"){
    quantumDot.MonteCarlo(task, whichMethod, alpha, beta, omega);
  }


  // Performing the parallelized MC sampling for hundred different alphas and
  // hundred different betas, with one set value for omega.
  // Creating one file for each alpha, containing hundred values for different betas
  else if (task == "Loop"){
    string filename;
    double beta0 = 0.20;
    double deltaBeta = 0.01;

    // Timing the parallelized code
    double start = omp_get_wtime();

    int n = 100;
    // Looping over the different values of alpha
    #pragma omp parallel for shared(n) private(filename)
    for (int i = 0; i < n; i++){
      double loopAlpha = alpha0 + i*deltaAlpha; // Updating alpha
      string alpha_, omega_;

      ostringstream streamObj1, streamObj2;
      streamObj1 << fixed << setprecision(2) << loopAlpha;
      streamObj2 << fixed << setprecision(2) << omega;
      alpha_ = streamObj1.str();
      omega_ = streamObj2.str();

      // Defining the filename
      filename = "../output/EnergyParallellized_";
      filename.append(alpha_).append("_").append(omega_).append(".dat");
      ofstream ofile;
      ofile.open(filename);

      // Looping over the different values of beta
      for (int j = 0; j < 100; j++){
        double loopBeta = beta0 + j*deltaBeta; // Updating beta

        quantumDot.MonteCarlo(task, whichMethod, loopAlpha, loopBeta, omega); // MC sampling
        quantumDot.WriteToFileParallel(ofile); // Writing values to file
      }
      quantumDot.CloseFile(ofile); // Closing file after use
    }

    // Finishing the timing
    double end = omp_get_wtime();
    double timeused = end-start;
    cout << "Parallelized code took " << timeused << " seconds to execute with number of MCCs = " << MCCs << "." << endl;
  }


  // Performing the MC sampling for hundred different values of the parameter omega
  else if ((task == "VirialwithoutInteraction") || (task == "VirialwithInteraction")){
    string filename;
    double omega;

    // Defining the filename
    filename = "../output/Energy_";
    filename.append(task).append(".dat");

    ofstream ofile;
    ofile.open(filename); // Opening the output file

    // Looping over the different values of omega
    for (int i = 0; i < 100; i++){
      omega = omega0 + i*deltaOmega; // Updating omega

      quantumDot.MonteCarlo(task, whichMethod, alpha, beta, omega); // MC sampling
      quantumDot.WriteToFileVirial(ofile); // Writing values to file
    }
    quantumDot.CloseFile(ofile); // Closing file after use
  }

  return 0;
}
