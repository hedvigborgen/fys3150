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

# To execute numerical solutions, and write the reults into texfile
./main.exe pow arg > 'textfile.txt'

#To compile relative error 
c++ -o error.exe error.cpp

#To execute relative error
./error.exe
```
Here; 
* pow is the power of ten;
	* which we assign as number of steps n 
* arg is the choise of algorithm; 
	* enter 1 for Gauss elimination for the general case
	* enter 2 for Gauss elimination for the special case
	* enter 3 for LU decomposition
	
To make plots of the results:
```
python3 readplot.py arg
```
Here;
* arg tells the program whether you want to make plots of the solution of the linear equation or its relative error; 
	* enter 1 to compare numerical solution with closed-form solution
	* enter 2 to compute the relative error for Gauss elimination for special case


## General info
This project is simple Lorem ipsum dolor generator.
	
	
## Setup
To run this project, install it locally using npm:

```
$ cd ../lorem
$ npm install
$ npm start
```
