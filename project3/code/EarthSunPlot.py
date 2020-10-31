import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import sys


# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


timestep = input('Enter number of time steps: ')
dt = input('Enter value for time step: ')

pos_files = ['euler', 'verlet']
energy_files = ['euler_system', 'verlet_system']
names = ['forward Euler', 'velocity Verlet']


for i in range(len(pos_files)):
    subprocess.call(['./main.exe', str(timestep), str(dt), str(i+1)])

def read_positions(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    n = int((len(lines))/2)
    x_Sun, y_Sun, z_Sun = np.zeros(n), np.zeros(n), np.zeros(n)
    x_Earth, y_Earth, z_Earth = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(0,len(lines)-1,2):
        line = lines[i]
        vals = line.split()
        x_Sun[i//2], y_Sun[i//2], z_Sun[i//2] = float(vals[1]), float(vals[2]), float(vals[3])
        line = lines[i+1]
        vals = line.split()
        x_Earth[i//2], y_Earth[i//2], z_Earth[i//2] = float(vals[1]), float(vals[2]), float(vals[3])
    infile.close()
    return x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n


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

for i, e in enumerate(pos_files):
    x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth, n = read_positions(f'../output/{e}.xyz')

    timestep = n
    fig, ax = plt.subplots()
    ax.plot(x_Sun, y_Sun, '*', color='#FFC800', label='Position of the Sun')
    ax.plot(x_Earth, y_Earth, color='#3498DB', label='Position of the Earth')
    plt.legend(fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_xlabel(r'x(t) [AU]', fontsize=15)
    ax.set_ylabel(r'y(t) [AU]', fontsize=15)
    ax.set_title(r'The Earth orbiting the Sun, '+'\n'+r'using {} method with n = {}'.format(names[i],timestep), fontsize=20)
    ax.axis('equal')
    fig.tight_layout()
    fig.savefig(f'../plots/positions_{e}{timestep}.pdf')



tot_E = np.zeros((2, timestep)) # Total energy as a function of time for both integration methods

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
    ax.set_title(r'Energies as a function of time, '+'\n'+r'using {} method with n = {}'.format(names[i],timestep), fontsize=20)
    fig.tight_layout()
    fig.savefig(f'../plots/energies_{e}{timestep}.pdf')

"""
# Computing relative change in total energy
rel_E_e = np.zeros(len(tot_E[0])-1)
rel_E_v = np.zeros(len(tot_E[1])-1)

for i in range(len(tot_E[0])-1):
    rel_E_e[i] = np.abs(tot_E[0, i+1] - tot_E[0, i])
for i in range(len(tot_E[1])-1):
    rel_E_v[i] = np.abs(tot_E[1, i+1] - tot_E[1, i])

"""

# Plotting absolute relative change in energy
fig, ax = plt.subplots()
ax.plot(t, tot_E[0], label='Forward Euler method')
ax.plot(t, tot_E[1], label='Velocity Verlet method')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r't [yr]', fontsize=15)
ax.set_ylabel(r'E$_{\text{tot}}$' + r'[$\text{M}_{\odot}\text{AU}^2/\text{yr}$]', fontsize=15)
ax.set_title(r'The total energy of the system as a function of time, '+'\n'+r' with n = {}'.format(timestep), fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/tot_energies_{timestep}.pdf')
