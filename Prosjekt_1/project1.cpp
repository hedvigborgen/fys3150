#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
#include "time.h"
#include <new>
#include <cstdio>
#include <cstring>
#include "lib.h"
using namespace std;


// Defining exact solution U(x) to the differential equation
double U(double x)
{
  return 1 - (1 - exp(-10))*x - exp(-10*x);
}


// Defining RHS of differential equation
double f(double x)
{
  return 100*exp(-10*x);
}


// Defining function for LU decomposition, finding inverse of matrix
void inverse(int n, double **matr)
{
int *index;
double d, *col, **y;

index = new int[n];
col  = new double[n];

 //for (int i = 0; i < n; i++) {
   //for (int j = 0; j< n; j++) {
     //cout << matr[i][j] << " ";
   //}
  // cout << endl;
 //}

y = (double **) matrix(n, n, sizeof(double));

ludcmp(matr, n, index, &d);   // LU decompose  a[][]

// find inverse of a[][] by columns
for(int j = 0; j < n; j++) {

// initialize right-side of linear equations
for(int i = 0; i < n; i++) col[i] = 0.0;
col[j] = 1.0;

lubksb(matr, n, index, col);

// save result in y[][]
for(int i = 0; i < n; i++) y[i][j] = col[i];
}   //j-loop over columns

// return the inverse matrix in a[][]
for(int i = 0; i < n; i++) {
  for(int j = 0; j < n; j++) matr[i][j] = y[i][j];
}

//for (int i = 0; i < n; i++) {
  //for (int j = 0; j< n; j++) {
    //cout << matr[i][j] << " ";
  //}
  //cout << endl;
//}

free_matrix((void **) y);     // release local memory
delete [] col;
delete []index;
}  // End: function inverse()






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
    clock_t start, finish;
    start = clock();

// Forward substitution, general case
  for (int i = 0; i < n-1; i++)
  {
    b_tilde[i+1] = b[i+1] - a[i]*c[i]/b_tilde[i];
    g_tilde[i+1] = g[i+1] - a[i]*g_tilde[i]/b_tilde[i];
  }


  double *v = new double[n];
  v[-1] = g_tilde[-1]/b_tilde[-1];

  // Backward substitution, general case
  for (int j = n-1; j > 0; j--)
  {
    v[j-1] = (g_tilde[j-1] - c[j-1]*v[j])/b_tilde[j-1];
  }

// Ending timer
  finish = clock();
  double time = (double (finish - start)/CLOCKS_PER_SEC);


// Output of x and y values
  cout << n << endl;
  for (int i = 0; i < n; i++)
  {
    cout << x[i] << " " << v[i] << endl;
  }

// Printing time for substitution
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

// Forward substitution, special case
for (int i = 1; i < n; i++)
{
  y[i] += y[i-1]*i/(i + 1);
}

y[n-1] /= (n + 1)/n;

// Backward substitution, special case
for (int j = n-2; j > -1; j--)
{
  y[j] = (y[j] + y[j+1])*(j+1)/(j+2);
}

// Ending timer
finish = clock();
double time = (double (finish - start)/CLOCKS_PER_SEC);


// Output of x and y values
cout << n << endl;
for (int i = 0; i < n; i++)
{
  //cout << x[i] << " " << y[i] << endl;
}

// Printing time for substitution
//cout << "Substitution took " << time << " seconds to execute in the special case." << endl;



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

//cout << eps_max << endl;



// Declaring matrix for LU decomposition
double **matr;

matr = (double **) matrix(n, n, sizeof(double));
for(int i = 0; i < n; i++)
{
for(int j = 0; j < n; j++)
{
  matr[i][j] = 0;
}
}

// Defining matrix A
matr[0][0] = 2;
for (int i = 1; i < n; i++)
{
  matr[i][i] = 2;
  matr[i-1][i] = -1;
  matr[i][i-1] = -1;
 }

// Finding inverse of matrix A
inverse(n, matr);

double *g = new double[n];
for (int i=0; i<n; i++)
{
  g[i] = pow(h,2)*f(x[i]);
}


// Matrix-vector multiplication A^{-1}*g = v
double *v = new double[n];
for (int i = 0; i < n; i++)
{
  for (int j = 0; j < n; j++)
  {
    v[i] += matr[i][j]*g[j];
  }
}

// Output of x and v values for LU decomposition
for (int i = 0; i < n; i++)
{
  cout << x[i] << " " << v[i] << endl;
}
delete [] x;
delete [] y;
delete [] u;
delete [] v;
delete [] g;

return 0;
}
}
