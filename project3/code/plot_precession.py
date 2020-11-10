import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining functions to read data files
def read_positions(filename, numberOfBodies):

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    numTimesteps = int((len(lines))/(numberOfBodies))

    positions = np.zeros((numTimesteps, numberOfBodies, 3))
    for i in range(numTimesteps):
        for j in range(numberOfBodies):
            line = lines[i*numberOfBodies + j]
            vals = line.split()
            positions[i, j] = [float(val) for val in vals[2:]]

    return positions


# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', '3', '1', '1', '1', '1', '1'])

precession_positions = read_positions('../output/verlet_positions.xyz', 2)


# Making the plot
fig, ax = plt.subplots()
ax.plot(precession_positions[:, 0, 0], precession_positions[:, 0, 1], '*', color='#FFC800', label='Sun')
ax.plot(precession_positions[:, 1, 0], precession_positions[:, 1, 1], color='#CD5C5C', label='Mercury')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'x(t) [AU]', fontsize=15)
ax.set_ylabel(r'y(t) [AU]', fontsize=15)
ax.set_title(f'The precession of Mercury for a century of its orbit', fontsize=20)
plt.legend(fontsize=15)
ax.axis('equal')
fig.tight_layout()
fig.savefig(f'../plots/precession.pdf')


arcsec = 206264.806 # 1 rad = 206264.806 arcsec

# Calculating the perihelion angle for precession with relativistic correction
theta_start = np.arctan(precession_positions[0, 1, 1]/precession_positions[0, 1, 0]) # Should be equal to zero for the initial position
theta_end = np.arctan(precession_positions[-1, 1, 1]/precession_positions[-1, 1, 0])

print(theta_end*arcsec)

# Calculating the perihelion angle for precession with a pure Newtonian force
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', '1', '24000', '0.001', '../input/precession_mercury.txt', '2', '2'])

newtonian_positions = read_positions('../output/verlet_positions.xyz', 2)
theta_newt_s = np.arctan(newtonian_positions[0, 1, 1]/newtonian_positions[0, 1, 0]) # Should be equal to zero for the initial position
theta_newt_e = np.arctan(newtonian_positions[-1, 1, 1]/newtonian_positions[-1, 1, 0])

print(theta_newt_e*arcsec)
