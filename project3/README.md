# Computational Physics, Project 3
## Simulation of the Solar system

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project, the numerical integration methods utilized are the forward Euler and velocity Verlet algorithms for solving differential equations governing the motions of celestial bodies in our Solar system.

## The project is created with:
* C++
   * iostream
   * math.h
   * vector
   * string
   * fstream
   * cmath
   * sstream
   * cstdlib
    
* Python version: 3
  * numpy	
  * matplotlib
  * subprocess
  * sys
  
* LaTeX

## How to run the code:
Open terminal window, these commands compile and execute the programs: 
```
# To compile the main script:
make compile

# To execute the same script:
./main.exe choice numTimesteps dt fname numberOfBodies method

```
Here; 
* choice determines the force of gravity;
    * Enter 1 for the pure Newtonian force
    * Enter 2 to test for different values of the parameter beta
    * Enter 3 to add relativistic correction to the Newtonian force
* numTimesteps determines the number of iterations;
* dt determines size of time step;
* fname determines the input file containing initial information;
* numberOfBodies determines the number of celestial bodies simulated;
* method determines the integration algorithm;
    * Enter 1 for forward Euler
    * Enter 2 for velocity Verlet
    
    
## To compile and execute the script;
```
# producing output files for ten celestial bodies by the velocity Verlet algorithm,
# with numTimesteps = 100,000, dt = 0.001:

make example

# producing output files with various values of beta by the velocity Verlet algorithm,
# with numTimesteps = 1000, dt = 0.001:

make compile
make test_beta

# producing output files for the precession of Mercury, numTimesteps = 1,000,000, dt = 0.0001:

make compile
make precession
```

### To make plots of positions of celestial bodies and mechanical energies for the solar system, with both forward Euler and velocity Verlet algorithm:
```
python3 plot_positions_energies.py numTimesteps dt numberOfBodies fname
```


### To make plots of positions of celestial bodies with force of gravity tested with various values of beta, with velocity Verlet algorithm:
```
python3 plot_test_forces.py
```


### To make a plot of the total angular momentum as a function of time for both elliptical and circular orbits with two-body system using velocity Verlet algorithm:
```
python3 plot_angmom.py numTimesteps dt
```

	
### To make a plot of the precession of Mercury's orbit, force of gravity with relativistic correction calulated with velocity Verlet method:
```
python3 plot_precession.py
```



## Example run 1: 
```
>> make compile                         				# Compiles the main script

>> ./main.exe 1 1000 0.001 ../input/two_bodies_elliptical.txt 2 2       # Produces output files for the positions of the
									celestial bodies, energy and angular momentum
									for the system, timing the integration

```


## Example run 2: Testing the force of gravity for various beta
```
>> python3 plot_test_forces.py                              

<< Enter number of time steps: 

>> 1000

<< Enter value for time step:

>> 0.001				# Compiles and executes the main script producing output files, makes plots
```

## Example run 3: Plotting the precession of Mercury's orbit
```
>> python3 plot_precession.py                                           

<< Enter number of time steps: 

>> 1000000

<< Enter dt:

>> 0.0001				# Compiles and executes the main script producing output files, makes plots

```

