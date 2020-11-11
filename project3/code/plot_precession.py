import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining function to read data files
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

    return positions, numTimesteps


# Defining function to find perihelion angle
def find_theta(positions, numTimesteps):
    strt = int(numTimesteps*0.99)
    dist = np.linalg.norm(positions[strt:, 1], axis=1)
    idx = np.where(dist == np.amin(dist))[0][0]
    peri_pos = positions[strt+idx, 1]

    theta_s = np.rad2deg(np.arctan2(positions[0, 1, 1], positions[0, 1, 0])) # Should be equal to zero for the initial position
    theta_e = np.rad2deg(np.arctan2(peri_pos[1], peri_pos[0]))

    return theta_s*3600, theta_e*3600


choice = '3'
n = input('Enter number of timesteps: ')
dt = input('Enter dt: ')

# Compiling and executing c++ script for gravitational force with relativistic correction
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', choice, n, dt, '1', '1', '1'])

# Calculating the perihelion angle for precession with relativistic correction
precession_positions, numTimesteps = read_positions('../output/verlet_positions.xyz', 2)
theta_s_gr, theta_e_gr = find_theta(precession_positions, numTimesteps)


# Compiling and executing c++ script for Newtonian force
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', '1', n, dt, '../input/precession_mercury.txt', '2', '2'])

# Calculating the perihelion angle for precession with a pure Newtonian force
newtonian_positions, numTimesteps = read_positions('../output/verlet_positions.xyz', 2)
theta_s_n, theta_e_n = find_theta(newtonian_positions, numTimesteps)

print(theta_e_n, theta_e_gr)


# Plotting the precession
fig, ax = plt.subplots()
ax.plot(precession_positions[:, 0, 0], precession_positions[:, 0, 1], '*', color='#FFC800', label='Sun')
ax.plot(precession_positions[:, 1, 0], precession_positions[:, 1, 1], color='#CD5C5C', label='Mercury')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'x(t) [AU]', fontsize=15)
ax.set_ylabel(r'y(t) [AU]', fontsize=15)
ax.set_title(f'The precession of Mercury for a century of its orbit', fontsize=20)
#plt.legend(fontsize=15)
ax.axis('equal')
fig.tight_layout()
fig.savefig(f'../plots/precession.pdf')
