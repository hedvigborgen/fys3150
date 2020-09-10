# Computational Physics, Project 1 
## Solving the one dimensional Poisson equation

Kamilla Ida Julie Sulebakk, Andrea Jensen Marthinussen and Hedvig Borgen Reiersrud

In this project we are solving a one dimensional Poisson equation. Here we will explore the numerical errors as we compare our numerical solution of the differential equation with an analytical solution. Using Dirichlet boundary conditions we will rewrite the differential equation as a set of linear equations, and we will thereafter use different linear algebra methods to solve it. The report is found above. 

## Project is created with:
* CPP
* Pyhton version: 3
	* Numpy 
	* Matplotlib

## How to run the code
Open terminal window and compile the code: 
```
#To compile numerical solutions 
c++ -o main.exe project1.cpp lib.cpp

# To execute numerical solutions, and write the results to file
./main.exe pow arg > 'textfile.txt'

#To compile relative error 
c++ -o error.exe error.cpp

#To execute relative error
./error.exe
```
Here; 
* pow is the power of ten;
	* which we assign as number of steps n 
* arg is the choice of algorithm; 
	* enter 1 for Gaussian elimination for the general case
	* enter 2 for Gaussian elimination for the special case
	* enter 3 for LU decomposition

	
To make plots of the results:
```
python3 readplot.py arg
```
Here;
* arg tells the program whether you want to make plots of the solution of the linear equation or its relative error; 
	* enter 1 to compare numerical solution with closed-form solution
	* enter 2 to compute the relative error for Gaussian elimination for special case
	
In the terminal window "Enter filename: " will be printed, in which case you would enter the name of the file containing the results you want to utilize. 

## Example run 1
```
c++ -o main.exe project1.cpp lib.cpp

./main.exe 3 2 > results_3_fast.txt   # runs the Gaussian elimination for the special case, with n = 10^(3)

python3 readplot.py 1                 # comparing the numerical solution with closed-form solution

>> Enter filename: 
<< results_3_fast.txt                 # produces a plot of the results in which you want to investigate
```

## Example run 2
```
c++ -o error.exe error.cpp

./error.exe			      

python3 readplot.py 2                 # plots the relative error for Gaussian elimination for special case
```
