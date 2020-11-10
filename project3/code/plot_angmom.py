import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining constant parameters
try:
    timestep, dt = sys.argv[1:]
    timestep = int(timestep)
    dt = float(dt)
except:
    timestep = int(input('Enter number of time steps: '))
    dt = float(input('Enter value for time step: '))

numberOfBodies = '2'
choice = '1'

# Defining function to read data file
def read_angmom(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    t, angmom = np.zeros(len(lines)),np.zeros(len(lines))

    for i in range(len(lines)):
        line = lines[i]
        vals = line.split()
        t[i], angmom[i] = float(vals[0]), float(vals[1])

    infile.close()
    return t, angmom

# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', choice, str(timestep), str(dt), '../input/two_bodies_circular.txt', numberOfBodies, '2'])

t, angmom_circular = read_angmom(f'../output/verlet_angmom.dat')

subprocess.call(['c++', '-o', 'main.exe', 'celestialbody.cpp', 'euler.cpp', 'main.cpp', 'mainfunc.cpp', 'solarsystem.cpp', 'vec3.cpp', 'velocityverlet.cpp', '--std=c++11'])
subprocess.call(['./main.exe', choice, str(timestep), str(dt), '../input/two_bodies_elliptical.txt', numberOfBodies, '2'])

t, angmom_elliptical = read_angmom(f'../output/verlet_angmom.dat')

# Plotting the angular momentum as a function of time, for both circular and elliptical orbits
fig, ax = plt.subplots()
ax.plot(t[1:], angmom_circular[1:], label='circular orbit')
ax.plot(t[1:], angmom_elliptical[1:], label='elliptical orbit')
ax.legend(fontsize=15)
ax.set_xlabel(r't [years]', fontsize=15)
ax.set_ylabel(r'L [M$_{\odot}\text{AU}^2$/yr]', fontsize=15)
ax.set_title('Angular momentum as a function of time', fontsize=20)
ax.set_ylim([0, 15])
fig.tight_layout()
fig.savefig(f'../plots/angmom_n{timestep}_dt{dt}.pdf')
