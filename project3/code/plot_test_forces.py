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
    n = int((len(lines) - 1)/2)
    x_Sun, y_Sun, z_Sun = np.zeros(n), np.zeros(n), np.zeros(n)
    x_Earth, y_Earth, z_Earth = np.zeros(n), np.zeros(n), np.zeros(n)
    beta = float(lines[0])
    for i in range(1,len(lines)-2,2):
        line = lines[i]
        vals = line.split()
        x_Sun[i//2], y_Sun[i//2], z_Sun[i//2] = float(vals[1]), float(vals[2]), float(vals[3])
        line = lines[i+1]
        vals = line.split()
        x_Earth[i//2], y_Earth[i//2], z_Earth[i//2] = float(vals[1]), float(vals[2]), float(vals[3])

    infile.close()
    return x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n, beta


# Defining arguments for c++ script
timestep = input('Enter number of time steps: ')
dt = input('Enter value for time step: ')
tot_time = timestep*dt
choice = '2'


# Data filenames
euler_files = ['euler_test_positions']
verlet_files = ['verlet_test_positions']


# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'main.cpp', 'celestialbody.cpp', 'solarsystem.cpp', 'euler.cpp', 'velocityverlet.cpp', 'vec3.cpp', '--std=c++11'])
methods = ['1', '2']
for i in range(len(methods)):
    subprocess.call(['./main.exe', str(timestep), str(dt), methods[i], choice])


#
for i, e in enumerate(euler_files):
    x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n, beta = read_positions(f'../output/{e}.xyz')

    fig, ax = plt.subplots()
    ax.plot(x_Sun, y_Sun, '*', color='#FFC800', label='Position of the Sun')
    ax.plot(x_Earth, y_Earth, color='#3498DB', label='Position of the Earth')
    plt.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'x(t) [AU]', fontsize=15)
    ax.set_ylabel(r'y(t) [AU]', fontsize=15)
    ax.set_title(r'The Earth orbiting the Sun, beta = {} '.format(beta) + '\n' + r'Forward Euler method for a time = {} years'.format(tot_time), fontsize=20)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/{e}_n{timestep}_dt{dt}.pdf')


#
for i, e in enumerate(verlet_files):
    x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n, beta = read_positions(f'../output/{e}.xyz')

    fig, ax = plt.subplots()
    ax.plot(x_Sun, y_Sun, '*', color='#FFC800', label='Position of the Sun')
    ax.plot(x_Earth, y_Earth, color='#3498DB', label='Position of the Earth')
    plt.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'x(t) [AU]', fontsize=15)
    ax.set_ylabel(r'y(t) [AU]', fontsize=15)
    ax.set_title(r'The Earth orbiting the Sun, beta = {} '.format(beta) + '\n' + r'Velocity Verlet method for a time = {} years'.format(tot_time), fontsize=20)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/{e}_n{timestep}_dt{dt}.pdf')
