#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
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
  int n = atof(argv[1]);
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

//Forward substitution:
for (int i = 1; i < n; i++)
{
  y[i] += y[i-1]*i/(i + 1); // Oppdaterer til f_tilde
}


y[n-1] /= (n + 1)/n;


// Backward substitution:
for (int j = n-2; j > -1; j--)
{
  y[j] = (y[j] + y[j+1])*(j)/(j+1); // Oppdaterer til u
}

cout << n << endl;
for (int i = 0; i < n; i++)
{
  cout << x[i] << " " << y[i] << endl;
}

delete [] x;
delete [] y;

return 0;

}
