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
    timestep, dt, numberOfBodies, file = sys.argv[1:]
    timestep = int(timestep)
    dt = float(dt)
    numberOfBodies = int(numberOfBodies)
except:
    timestep = int(input('Enter number of time steps: '))
    dt = float(input('Enter value for time step: '))
    numberOfBodies = int(input('Enter number of bodies: '))
    if numberOfBodies == 2:
        orbit = int(input("Enter 1 for circular orbit, or 2 for elliptical orbit: "))
        if orbit == 1:
            file = "../input/two_bodies_circular.txt"
        elif orbit == 2:
            file = "../input/two_bodies_elliptical.txt"
    elif numberOfBodies == 3:
        file = "../input/three_bodies.txt"
    elif numberOfBodies == 10:
        file = "../input/ten_bodies.txt"


# method = input('Enter 1 for Euler, 2 for Verlet, or 3 for both: ')
choice = '1'

# Data filenames
pos_input = ['euler_positions', 'verlet_positions']
energy_input = ['euler_energies', 'verlet_energies']


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
subprocess.call(['c++','-o','main.exe','$(wildcard *.cpp)','--std=c++11'])
for i in range(len(pos_input)):
    subprocess.call(['./main.exe', choice, str(timestep), str(dt), str(file), str(numberOfBodies), str(i+1)])
#
# if method == '1':
#     subprocess.call(['./main.exe', str(timestep), str(dt), str(1), str(file), str(numberOfBodies), choice])
#
# elif method == '2':
#     subprocess.call(['./main.exe', str(timestep), str(dt), str(2), str(file), str(numberOfBodies), choice])
#
# elif method == '3':
#     for i in range(len(pos_input)):
#         subprocess.call(['./main.exe', str(timestep), str(dt), str(i+1), str(file), str(numberOfBodies), choice])


methods = ['forward Euler', 'velocity Verlet']
colors = ['#3498DB', '#33D7FF', '#9EC1CF', '#9EE09E', '#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']#, '#EE452A']

n_str = f'$n = {timestep:.1e}'.replace('e+0', r'\cdot 10^{').replace('e-0', r'\cdot 10^{-') + r'}$'
dt_str = f'$dt = {dt:.1e}'.replace('e+0', r'\cdot 10^{').replace('e-0', r'\cdot 10^{-') + r'}$'
title = '\n with ' + n_str + ' and ' + dt_str
outfile = f'{file}'.replace('../input/','').replace('.txt', '')


for i, pos_file in enumerate(pos_input):
    names, time, positions, numTimesteps = read_positions(f'../output/{pos_file}.xyz', numberOfBodies)
    dt = time[-1]/(numTimesteps-1)
    fig, ax = plt.subplots()
    ax.plot(positions[:, 0, 0], positions[:, 0, 1], '*', color='#FFC800', label='Sun')

    for j,f in enumerate(names[1:]):
        ax.plot(positions[:, j+1, 0], positions[:, j+1, 1], color=colors[j], label=f'{f}')
        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.set_xlabel(r'x(t) [AU]', fontsize=15)
        ax.set_ylabel(r'y(t) [AU]', fontsize=15)

    # black = plt.gca()
    # black.set_facecolor('black')
    # l = plt.legend(fontsize=15)
    # for text in l.get_texts():
    #     text.set_color("#89B5C8")
    ax.set_title(f'The Earth orbiting the Sun using {methods[i]} method, ' + title, fontsize=20)
    plt.legend(fontsize=15)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/{pos_file}_n{timestep}_dt{dt}_{outfile}.pdf')



tot_E = np.zeros((2, int(timestep))) # Total energy as a function of time for both integration methods

for i, energy_file in enumerate(energy_input):
    t, pot_E, kin_E, tot_E[i] = read_energies(f'../output/{energy_file}.dat')

    # Plotting the energy as a function of time
    fig, ax = plt.subplots()
    ax.plot(t[1:], pot_E[1:], label='Potential energy')
    ax.plot(t[1:], kin_E[1:], label='Kinetic energy')
    ax.plot(t[1:], tot_E[i, 1:], label='Total energy')
    ax.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r't [years]', fontsize=15)
    ax.set_ylabel(r'E [M$_{\odot}\text{AU}^2$/yr]', fontsize=15)
    ax.set_title(f'Energies as a function of time using {methods[i]} method, ' + title, fontsize=20)
    fig.tight_layout()
    fig.savefig(f'../plots/{energy_file}_n{timestep}_dt{dt}_{outfile}.pdf')


# Comparing absolute change in energy between Euler and Verlet method
fig, ax = plt.subplots()
ax.plot(t[1:], tot_E[0, 1:], label='Forward Euler method')
ax.plot(t[1:], tot_E[1, 1:], label='Velocity Verlet method')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r't [yr]', fontsize=15)
ax.set_ylabel(r'E$_{\text{tot}}$' + r'[$\text{M}_{\odot}\text{AU}^2/\text{yr}$]', fontsize=15)
ax.set_title(f'The total energy of the system as a function of time, ' + title, fontsize=20)
ax.set_ylim([-0.00025, 0.00005])
fig.tight_layout()
fig.savefig(f'../plots/compare_energies_n{timestep}_dt{dt}_{outfile}.pdf')
