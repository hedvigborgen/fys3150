import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining constant parameters
timestep = input('Enter number of time steps: ')
dt = input('Enter value for time step: ')
file = input('Enter input filename: ') #'../files/ten_bodies.txt'
numberOfBodies = int(input('Enter number of bodies: '))
choice = '1'

# Data filenames
pos_files = ['euler_positions', 'verlet_positions']
energy_files = ['euler_energies', 'verlet_energies']


# Defining functions to read data files
def read_positions(filename, numberOfBodies):

    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()

    numTimesteps = int((len(lines))/(numberOfBodies))

    positions = np.zeros((numTimesteps, numberOfBodies, 3))
    time = np.zeros(numTimesteps)
    names = []
    for i in range(numTimesteps):
        for j in range(numberOfBodies):
            line = lines[i*numberOfBodies + j]
            vals = line.split()
            if len(names) < numberOfBodies:
                names.append(vals[0])
            if j == 0:
                time[i] = vals[1]
            positions[i, j] = [float(val) for val in vals[2:]]

    return names, time, positions, numTimesteps



def read_energies(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    t, pot_E, kin_E, tot_E = np.zeros(len(lines)),np.zeros(len(lines)), np.zeros(len(lines)), np.zeros(len(lines))

    for i in range(len(lines)):
        line = lines[i]
        vals = line.split()
        t[i], pot_E[i], kin_E[i], tot_E[i] = float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])

    infile.close()
    return t, pot_E, kin_E, tot_E


# Compiling and executing c++ script
subprocess.call(['c++', '-o', 'main.exe', 'main.cpp', 'celestialbody.cpp', 'solarsystem.cpp', 'euler.cpp', 'velocityverlet.cpp', 'vec3.cpp', '--std=c++11'])
for i in range(len(pos_files)):
    subprocess.call(['./main.exe', str(timestep), str(dt), str(file), str(numberOfBodies), str(i+1), choice])


methods = ['forward Euler', 'velocity Verlet']
colors = ['#3498DB', '#FFC800', '#9EC1CF', '#9EE09E', '#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#CC99C9', '#EE452A']

for i, e in enumerate(pos_files):
    names, time, positions, numTimesteps = read_positions(f'../output/{e}.xyz', numberOfBodies)
    dt = time[-1]/(numTimesteps-1)
    fig, ax = plt.subplots()
    ax.plot(positions[:, 0, 0], positions[:, 0, 1], '*', color='#FFC800', label='Sun')

    for j,f in enumerate(names[1:]):
        ax.plot(positions[:, j+1, 0], positions[:, j+1, 1], color=colors[j], label=f'{f}')
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel(r'x(t) [AU]', fontsize=15)
        ax.set_ylabel(r'y(t) [AU]', fontsize=15)

    ax.set_title(r'The Earth orbiting the Sun, '+'\n'+r'using {} method with dt = {}'.format(methods[i],dt), fontsize=20)
    plt.legend(fontsize=15)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/positions_{e}_n{timestep}_dt{dt}_{numberOfBodies}.pdf')



tot_E = np.zeros((2, int(timestep))) # Total energy as a function of time for both integration methods

for i, e in enumerate(energy_files):
    t, pot_E, kin_E, tot_E[i] = read_energies(f'../output/{e}.dat')

    #Plotting the energy as a function of time
    fig, ax = plt.subplots()
    ax.plot(t, pot_E, label='Potential energy')
    ax.plot(t, kin_E, label='Kinetic energy')
    ax.plot(t, tot_E[i], label='Total energy')
    ax.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r't [years]', fontsize=15)
    ax.set_ylabel(r'E [M$_{\odot}\text{AU}^2$/yr]', fontsize=15)
    ax.set_title(r'Energies as a function of time, '+'\n'+r'using {} method with n = {}'.format(methods[i],timestep), fontsize=20)
    fig.tight_layout()
    fig.savefig(f'../plots/energies_{e}_n{timestep}_dt{dt}_{numberOfBodies}.pdf')


# Comparing absolute change in energy between Euler and Verlet method
fig, ax = plt.subplots()
ax.plot(t, tot_E[0], label='Forward Euler method')
ax.plot(t, tot_E[1], label='Velocity Verlet method')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r't [yr]', fontsize=15)
ax.set_ylabel(r'E$_{\text{tot}}$' + r'[$\text{M}_{\odot}\text{AU}^2/\text{yr}$]', fontsize=15)
ax.set_title(r'The total energy of the system as a function of time', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/compare_energies_n{timestep}_dt{dt}_{numberOfBodies}.pdf')
