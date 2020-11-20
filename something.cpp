for ( j=1; j < N; j++) {
energy += spin[j]*spin[j+1];
}

jm=N;
for ( j=1; j <=N ; j++) {
            energy += spin[j]*spin[jm];
jm = j ; }


#include <omp.h>


int main(int argc, char const *argv[]) {

  #pragma omp parallel for
  for (int i = 0; i < num_temps; i++){
    temp = T0 + i*dT;
    Ising2D my_solver(double temp, int L);

  }



  return 0;
}