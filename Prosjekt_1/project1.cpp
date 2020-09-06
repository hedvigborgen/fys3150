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

int arg = atoi(argv[2]);

// THE GENERAL CASE
if (arg == 0)
{
// Defining matrix elements for the general case
  double *a = new double[n-1];
  double *b = new double[n];
  double *b_tilde = new double[n];
  double *c = new double[n-1];

  double *g = new double[n];
  double *g_tilde = new double[n];

for (int i = 0; i < n-1; i++)
{
  a[i] = -1;
}

for (int i = 0; i < n-1; i++)
{
  c[i] = -1;
}

for (int i = 0; i < n; i++)
  {
    b[i] = 2;
  }

for (int i = 0; i < n; i++)
  {
    g[i] = pow(h, 2)*f(x[i]);
  }

  b_tilde[0] = b[0];
  g_tilde[0] = g[0];


// Timing the Gauss elimination for the general case
    clock_t start, finish;  //  declare start and final time
    start = clock();

// Forward substitution, general case:
  for (int i = 0; i < n-1; i++)
  {
    b_tilde[i+1] = b[i+1] - a[i]*c[i]/b_tilde[i];
    g_tilde[i+1] = g[i+1] - a[i]*g_tilde[i]/b_tilde[i];
  }


  double *v = new double[n];
  v[-1] = g_tilde[-1]/b_tilde[-1];

  // Backward substitution, general case:
  for (int j = n-1; j > 0; j--)
  {
    v[j-1] = (g_tilde[j-1] - c[j-1]*v[j])/b_tilde[j-1];
  }

  finish = clock();
  double time = (double (finish - start)/CLOCKS_PER_SEC);


  cout << n << endl;
  for (int i = 0; i < n; i++)
  {
    cout << x[i] << " " << v[i] << endl;
  }
  cout << "Substitution took " << time << " seconds to execute in the general case." << endl;
  delete [] a;
  delete [] b;
  delete [] b_tilde;
  delete [] c;
  delete [] g;
  delete [] g_tilde;
  delete [] v;
}


// THE SPECIAL CASE
else if (arg == 1)
{
// Defining matrix elements for the special case
  int a = -1;
  int b = 2;
  int c = -1;

  for (int i=0; i<n; i++)
  {
    y[i] = pow(h,2)*f(x[i]);
  }

// Timing the Gauss elimination for the special case
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
cout << "Substitution took " << time << " seconds to execute in the special case." << endl;


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

delete [] u;
}


delete [] x;
delete [] y;


return 0;
}
