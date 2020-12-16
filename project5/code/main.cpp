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

  dimension = 3;
  numberofParticles = 2;
  task = arguments[1];


  if (task == "MCCs"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    maxVariations = atol(arguments[4]);
    alpha0 = atof(arguments[5]);
    deltaAlpha = atof(arguments[6]);

    equilibrationTime = 0;
    write = "each time";
    step = 1.0;
    beta = 1.0;
    omega = 1.0;
    charge = 1.0;
  }

  else if (task == "StepLength"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    maxVariations = atol(arguments[4]);
    alpha0 = atof(arguments[5]);
    deltaAlpha = atof(arguments[6]);
    step = atof(arguments[7]);

    equilibrationTime = 100000;
    write = "at the end";
    beta = 1.0;
    omega = 1.0;
    charge = 1.0;
  }

  else if (task == "Parameters"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    whichMethod = atoi(arguments[3]);
    alpha = atof(arguments[4]);
    beta = atof(arguments[5]);
    omega = atof(arguments[6]);

    maxVariations = 1;
    equilibrationTime = 100000;
    charge = 1.0;

    if ((whichMethod != 0) && (whichMethod != 1) && (whichMethod != 2)){
      cout << "Unacceptable choice of method!" << endl; exit(1);
    }
  }

  else if (task == "Loop"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega = atof(arguments[3]);

    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 2;
    alpha0 = 0.60;
    deltaAlpha = 0.01;
    charge = 1.0;
  }

  else if (task == "VirialwithoutInteraction"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega0 = atof(arguments[3]);
    deltaOmega = atof(arguments[4]);

    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 0;
    alpha = 0.995;
    beta = 0.280;
    charge = 1.0;
  }

  else if (task == "VirialwithInteraction"){
    // Defining input arguments
    MCCs = atol(arguments[2]);
    omega0 = atof(arguments[3]);
    deltaOmega = atof(arguments[4]);

    maxVariations = 1;
    equilibrationTime = 100000;
    whichMethod = 2;
    alpha = 0.995;
    beta = 0.280;
    charge = 1.0;
  }



  else if ((task != "MCCs") && (task != "StepLength") && (task != "Parameters")
  && (task != "Loop") && (task != "VirialwithoutInteraction") && (task != "VirialwithInteraction")){
    cout << "Unacceptable choice of task!" << endl; exit(1);
  }


  // Initializing the system with input arguments
  QuantumDot quantumDot(dimension, numberofParticles, charge, equilibrationTime, MCCs);


  if (maxVariations != 1){
    // Performing the MC sampling
    quantumDot.MonteCarlo(whichMethod, write, maxVariations, alpha0, deltaAlpha, beta, omega, step);
  }


  else if (task == "Parameters"){
    // Performing the MC sampling
    quantumDot.MonteCarlo(task, whichMethod, alpha, beta, omega);
  }


  else if (task == "Loop"){
    string filename;
    double beta0 = 0.20;
    double deltaBeta = 0.01;

    // Timing the parallelized code
    double start = omp_get_wtime();

    // Performing the MC sampling
    #pragma omp parallel for private(filename)
    for (int i = 0; i < 100; i++){
      double loopAlpha = alpha0 + i*deltaAlpha;
      string alpha_, omega_;

      ostringstream streamObj1, streamObj2;
      streamObj1 << fixed << setprecision(2) << loopAlpha;
      streamObj2 << fixed << setprecision(2) << omega;
      alpha_ = streamObj1.str();
      omega_ = streamObj2.str();

      filename = "../output/EnergyParallellized_";
      filename.append(alpha_).append("_").append(omega_).append(".dat");
      ofstream ofile;
      ofile.open(filename);

      for (int j = 0; j < 100; j++){
        double loopBeta = beta0 + j*deltaBeta;

        quantumDot.MonteCarlo(task, whichMethod, loopAlpha, loopBeta, omega);
        quantumDot.WriteToFileParallel(ofile);
      }
      quantumDot.CloseFile(ofile);
    }

    // Finishing the timing
    double end = omp_get_wtime();
    double timeused = end-start;
    cout << "Parallelized code took " << timeused << " seconds to execute with number of MCCs = " << MCCs << "." << endl;
  }


  else if ((task == "VirialwithoutInteraction") || (task == "VirialwithInteraction")){
    string filename;
    double omega;

    filename = "../output/Energy_";
    filename.append(task).append(".dat");

    ofstream ofile;
    ofile.open(filename);

    // Performing the MC sampling
    for (int i = 0; i < 100; i++){
      omega = omega0 + i*deltaOmega;

      quantumDot.MonteCarlo(task, whichMethod, alpha, beta, omega);
      quantumDot.WriteToFileVirial(ofile);
    }
    quantumDot.CloseFile(ofile);
  }

  return 0;
}
