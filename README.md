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
./main.exe n arg > 'textfile.txt' # where n is the size of the number of , arg is 

#To compile relative error 
c++ -o error.exe error.cpp

#To execute relative error
./error.exe
```
To produce plots of numerical soultions, for a specific n, run: 
```
python3 readplot.py 
```


$ cd ../lorem
$ npm install
$ npm start



* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)

## General info
This project is simple Lorem ipsum dolor generator.
	
	
## Setup
To run this project, install it locally using npm:

```
$ cd ../lorem
$ npm install
$ npm start
```
