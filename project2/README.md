# Computational Physics, Project 2 
## Jacobi's method for solving differential equations

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project we use Jacobi’s method to solve three second order differential equations. The first equation is a classical problem, known as the buckling beam problem. The latter two are Schrödinger equations describing systems affected by an harmonic oscillator potential, with one and two electrons respectively.

## The project is created with:
* C++
  * iostream
  * math.h
  * armadillo
* Python version: 3
	* Numpy 
	* Matplotlib
  * Subprocess
* LaTeX

## How to run the code:
Open terminal window, these commands compile and execute the programs: 
```
#To compile the main script:
make compile

# To execute the same script:
./main.exe n rho_max omega_r method

#To compile and execute the script containing unit tests:
make all_test

```
Here; 
* n determines the size of the nxn-matrix A;
* rho_max is the dimensionless maximum value of the position relative to a reference point;
* omega_r is the parameter reflecting the strength of an oscillating potential;
* method determines your output; 
	* enter count to print number of iterations needed for diagonalization of matrix A
  * enter plotdiff to print numerical deviation from the smallest analytic eigenvalue of A
  * enter ploteig to print initial eigenvector for smallest eigenvalue of A
  * enter eigvals to print initial eigenvalues of matrix A
	* enter eigvecs to print eigenvectors of matrix A computed different ways

	
To plot the difference between the numerical and analytical values for the smallest eigenvalue of matrix A as a function of rho_max, for different sizes nxn of A:
```
python3 rho_vs_diff.py
```

To make plots of the analytical and the numerical solution u(rho) as a function of the dimensionless position on the beam $\rho$:
```
python3 eigenvector_plot.py
```

To make a plot of the probability |psi(rho)|^2 as a function of the dimensionless distance rho between two electrons for different omega_r:
```
python3 eigenvector_QM_plot.py
```
	


## Example run 1: Main script
```
>> make compile                                             # Compiles the main script

>> ./main.exe 100 5 0 count                                 # Runs main script with n = 100, rho_max = 5,
                                                              omega_r = 0 for method count

<< Number of iterations needed for diagonalization: 16475   # Prints number of iterations needed for
                                                              diagonalization of matrix A
```

## Example run 2: Test script
```
>> make all_test                                            # Compiles and executes the script containing unit tests

<<  Correct maximum value found.
Eigenvalues are preserved.                	            # Output is the result of the unit tests
```


## Example run 3: Ploting 
```
>> python3 eigenvector_QM_plot.py                           # Compiles the main script

<< resultsQM.pdf                                            # Produces plot of probability as a function of rho

```
