import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining arguments for c++ script
choice = '2'
timesteps = int(input('Enter number of time steps: '))
dt = float(input('Enter value for time step: '))
tot_time = timesteps*dt
filename = "../input/two_bodies_elliptical.txt"
num_bodies = '2'
method = '2'


# Defining function to read data files
def read_positions(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    n = timesteps
    time = np.linspace(0, tot_time, n)
    x_Sun, y_Sun, z_Sun = np.zeros((6,n)), np.zeros((6,n)), np.zeros((6,n))
    x_Earth, y_Earth, z_Earth = np.zeros((6,n)), np.zeros((6,n)), np.zeros((6,n))
    beta = np.zeros(6)

    k = 0
    for i in range(6):
        beta[i] = float(lines[k])

        k += 1
        for j in range(0, n):
            line = lines[k]
            vals = line.split()
            x_Sun[i, j], y_Sun[i, j], z_Sun[i, j] = float(vals[1]), float(vals[2]), float(vals[3])
            line = lines[k+1]
            vals = line.split()
            x_Earth[i, j], y_Earth[i, j], z_Earth[i, j] = float(vals[1]), float(vals[2]), float(vals[3])
            k += 2

    infile.close()
    return x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, beta - 1, time


# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', choice, str(timesteps), str(dt), filename, num_bodies, method])


# Plotting
x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, beta, time = read_positions('../output/verlet_test_positions.xyz')
for j in range(len(beta)):
    fig, ax = plt.subplots()
    ax.plot(y_Sun[j], z_Sun[j], '*', color='#FFC800', label='Position of the Sun')
    ax.plot(y_Earth[j], z_Earth[j], color='#3498DB', label='Position of the Earth')
    plt.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'y(t) [AU]', fontsize=15)
    ax.set_ylabel(r'z(t) [AU]', fontsize=15)
    ax.set_title(r'The Earth orbiting the Sun, beta = {} '.format(beta[j]) + '\n' + r'velocity Verlet method for a time = {} years'.format(tot_time), fontsize=20)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/plot_beta{beta[j]}_n{timesteps}_dt{dt}.pdf')
