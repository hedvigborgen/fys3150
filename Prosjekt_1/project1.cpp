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

int main ()
{
  int n = 100;
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
for (int i = 0; i < n-1; i++)
{
  y[i+1] += y[i]*(i + 1)/(i + 2); // Oppdaterer til f_tilde
}


y[n-1] /= (n + 1)/n;


// Backward substitution:
for (int j = n-1; j > 0; j--)
{
  y[j-1] += y[j]*(j+1)/(j+2); // Oppdaterer til u
}

for (int i = 0; i < n; i++)
{
  cout << x[i] << " " << y[i] << endl;
}

delete [] x;
delete [] y;

return 0;

}
