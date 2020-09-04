#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
#include "time.h"
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

int main (int argc, char** argv)
{
  int n = pow(10, atof(argv[1]));
  double *x = new double[n];
  double *y = new double[n];
  double h = 1/(n-1.0);

  for (int i=0; i<n; i++)
  {
    x[i] = i*h;
  }


  int a = -1;
  int b = 2;
  int c = -1;

  for (int i=0; i<n; i++)
  {
    y[i] = pow(h,2)*f(x[i]);
  }


// Timing the Gauss elimination for special case
  clock_t start, finish;  //  declare start and final time
  start = clock();

// Forward substitution, special case:
for (int i = 1; i < n; i++)
{
  y[i] += y[i-1]*i/(i + 1); // Oppdaterer til f_tilde
}


y[n-1] /= (n + 1)/n;


// Backward substitution, special case:
for (int j = n-2; j > -1; j--)
{
  y[j] = (y[j] + y[j+1])*(j)/(j+1); // Oppdaterer til u
}

finish = clock();
double time = (double (finish - start)/CLOCKS_PER_SEC);


cout << n << endl;
for (int i = 0; i < n; i++)
{
  cout << x[i] << " " << y[i] << endl;
}
cout << "Special substitution took " << time << " seconds to execute." << endl;

// Computing the relative error
double *u = new double[n];

for (int i=0; i<n; i++)
{
  u[i] = U(x[i]);
}

double eps_max, eps;
eps_max = 0;
for(int i=0; i<n; i++)
  {
    eps = log(abs((y[i]-u[i])/u[i]));
  }
  if (eps > eps_max)
    {
      eps_max = eps;
    }

cout << eps_max << endl;


delete [] x;
delete [] y;
delete [] u;

return 0;

}
