#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
#include "time.h"
#include <vector>
using namespace std;

// Defining exact solution:
double U(double x)
{
  return 1 - (1 - exp(-10))*x - exp(-10*x);
}
// Defining RHS:
double f(double x)
{

  return 100*exp(-10*x);

}

int main ()
{

  int j = 7;
  int *n = new int[j];
  for (int i=1; i<j+1; i++){
    n[i-1] = pow(10, i);
  }

  for (int k=0; k<j; k++){


    double *x = new double[n[k]];
    double *y = new double[n[k]];
    double h = 1/(n[k]-1.0);

    for (int i=0; i<n[k]; i++)
    {
      x[i] = i*h;
    }


    int a = -1;
    int b = 2;
    int c = -1;

    for (int i=0; i<n[k]; i++)
    {
      y[i] = pow(h,2)*f(x[i]);
    }




  // Forward substitution, special case:
  for (int i = 1; i < n[k]; i++)
  {
    y[i] += y[i-1]*i/(i + 1); // Oppdaterer til f_tilde
  }


  y[n[k]-1] /= (n[k] + 1)/n[k];


  // Backward substitution, special case:
  for (int j = n[k]-2; j > -1; j--)
  {
    y[j] = (y[j] + y[j+1])*(j)/(j+1); // Oppdaterer til u
  }



  // Computing the relative error
  double *u = new double[n[k]];

  for (int i=0; i<n[k]; i++)
  {
    u[i] = U(x[i]);
  }

  double eps_max, eps;
  eps_max = 0;
  for(int i=0; i<n[k]; i++)
    {
      eps = log(abs((y[i]-u[i])/u[i]));
    }
    if (eps > eps_max)
      {
        eps_max = eps;
      }

  cout << eps_max << endl;



  }


return 0;

}
