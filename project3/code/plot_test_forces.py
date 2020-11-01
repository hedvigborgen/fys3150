import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Defining function to read data files
def read_positions(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    n = timesteps
    x_Sun, y_Sun, z_Sun = np.zeros((6,n)), np.zeros((6,n)), np.zeros((6,n))
    x_Earth, y_Earth, z_Earth = np.zeros((6,n)), np.zeros((6,n)), np.zeros((6,n))
    beta = np.zeros(6)

    k = 0
    for i in range(6):
        beta[i] = float(lines[k])
        k += 1
        for j in range(0, 1000):
            line = lines[k]
            print(line)
            vals = line.split()
            x_Sun[i, j], y_Sun[i, j], z_Sun[i, j] = float(vals[1]), float(vals[2]), float(vals[3])
            line = lines[k+1]
            vals = line.split()
            x_Earth[i, j], y_Earth[i, j], z_Earth[i, j] = float(vals[1]), float(vals[2]), float(vals[3])
            k += 2

    infile.close()
    return x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, beta


# Defining arguments for c++ script
timesteps = int(input('Enter number of time steps: '))
dt = float(input('Enter value for time step: '))
tot_time = timesteps*dt
choice = '2'


# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'main.cpp', 'celestialbody.cpp', 'solarsystem.cpp', 'euler.cpp', 'velocityverlet.cpp', 'vec3.cpp', '--std=c++11'])
methods = ['1', '2']
for method in methods:
    subprocess.call(['./main.exe', str(timesteps), str(dt), method, choice])


# Lists of data filenames and method names
filenames = ['euler_test_positions', 'verlet_test_positions']
method_names = ['Forward Euler', 'Velocity Verlet']


for file in filenames:
    x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n, beta = read_positions(f'../output/{file}.xyz')
    for j in range(len(beta)):
        fig, ax = plt.subplots()
        ax.plot(x_Sun[j], y_Sun[j], '*', color='#FFC800', label='Position of the Sun')
        ax.plot(x_Earth[j], y_Earth[j], color='#3498DB', label='Position of the Earth')
        plt.legend(fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel(r'x(t) [AU]', fontsize=15)
        ax.set_ylabel(r'y(t) [AU]', fontsize=15)
        ax.set_title(r'The Earth orbiting the Sun, beta = {} '.format(beta[j]) + '\n' + r'{} method for a time = {} years'.format(method_names[filenames.index(file)], tot_time), fontsize=20)
        ax.axis('equal')
        fig.tight_layout()
        fig.savefig(f'../plots/plot_{file}_beta{beta[j]}_n{timesteps}_dt{dt}.pdf')
