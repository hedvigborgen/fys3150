#include "IsingModel.hpp"
#include <omp.h>


int main(int numArg, char *arguments[]) {
  string part, filename;
  int MCCs, L, whichMatrix;
  double T;

  part = arguments[1];
  MCCs = atoi(arguments[2]);
  L = atoi(arguments[3]);


  // Running the script for part 4c
  if (part == "c"){
    T = atof(arguments[4]);
    whichMatrix = atoi(arguments[5]);

    int ChoiceofSampling = 1;
    filename = "../output/part4c_";
    filename.append(arguments[4]).append(".dat");

    // Initializing the system
    IsingModel isingModel(L, whichMatrix, T);

    // Calculating observables for the spin matrix and writing to file
    isingModel.MetropolisCycle(MCCs, filename, whichMatrix, ChoiceofSampling);
  }


  // Running the script for part 4d
  if (part == "d"){
    T = atof(arguments[4]);
    whichMatrix = atoi(arguments[5]);
    int ChoiceofSampling = 1;

    if (whichMatrix == 1){
      filename = "../output/part4d_ordered_";
      filename.append(arguments[4]).append(".dat");
    }

    else if (whichMatrix == 2){
      filename = "../output/part4d_random_";
      filename.append(arguments[4]).append(".dat");
    }

    // Initializing the system
    IsingModel isingModel(L, whichMatrix, T);

    // Calculating observables for the system and writing to file
    isingModel.MetropolisCycle(MCCs, filename, whichMatrix, ChoiceofSampling);
  }


  // Running the script for part 4e
  if (part == "e"){
    T = atof(arguments[4]);
    int ChoiceofSampling = 2;
    whichMatrix = 2;
    filename = "../output/part4e_random_";
    filename.append(arguments[4]).append(".dat");

    // Initializing the system
    IsingModel isingModel(L, whichMatrix, T);

    // Calculating observables for the system and writing to file
    isingModel.MetropolisCycle(MCCs, filename, whichMatrix, ChoiceofSampling);
  }


  // Running the script for part 4f
  else if (part == "f"){
    whichMatrix = atoi(arguments[4]);

    double T_start = 2.15;
    double T_end = 2.45;
    int Ntemps = 20;
    double h = (T_end - T_start)/(Ntemps - 1);


    // Timing the parallelized code
    double start = omp_get_wtime();

    #pragma omp parallel for private(filename)
    for (int i = 0; i < Ntemps; i++){
      T = T_start + i*h;
      filename = "../output/part4f_";
      filename.append(to_string(i)).append("_L_").append(arguments[3]).append(".dat");

      // Initializing the system
      IsingModel isingModel(L, whichMatrix, T);

      // Calculating observables for the system and writing to file
      isingModel.MetropolisCycle(MCCs, filename);
    }

    // Finishing the timing
    double end = omp_get_wtime();
    double timeused = end-start;
    cout << "Parallelized code took " << timeused << " seconds to execute with L = " << L << "." << endl;
  }


else if (part == "unparallelized_f"){
  whichMatrix = atoi(arguments[4]);

  double T_start = 2.15;
  double T_end = 2.45;
  int Ntemps = 20;
  double h = (T_end - T_start)/(Ntemps - 1);

  // Timing the test
  clock_t start, finish;
  start = clock();

  for (int i = 0; i < Ntemps; i++){
    T = T_start + i*h;
    filename = "../output/test_";
    filename.append(to_string(i)).append("_L_").append(arguments[3]).append(".dat");

    // Initializing the system
    IsingModel isingModel(L, whichMatrix, T);

    // Calculating observables for the system and writing to file
    isingModel.MetropolisCycle(MCCs, filename);
  }

  // Finishing the timing
  finish = clock();
  double time = (double (finish - start)/CLOCKS_PER_SEC);
  cout << "Unparallelized code took " << time << " seconds to execute with L = " << L << "." << endl;

}

  return 0;
}
