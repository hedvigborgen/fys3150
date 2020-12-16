import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Defining input arguments
tasks = ['VirialwithoutInteraction', 'VirialwithInteraction']
MCCs = 1_000_000
omega0 = 0.01
deltaOmega = 0.01
size = 100


omega = np.linspace(omega0, deltaOmega*size, size)
# Arrays for storing the energies
kineticEnergy = np.zeros((len(tasks), size))
potentialEnergy = np.zeros((len(tasks), size))


# Compiling and executing the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
for i, task in enumerate(tasks):
    # Reading the output files
    subprocess.call(['./main.exe', task, f'{MCCs}', f'{omega0}', f'{deltaOmega}'])
    filename = f'../output/Energy_{task}.dat'
    infile = open(filename, 'r')
    lines = infile.readlines()

    # Storing the important values
    for j, line in enumerate(lines):
        vals = line.split()
        kineticEnergy[i,j] = float(vals[1])
        potentialEnergy[i,j] = float(vals[2])

    infile.close()


# Calculating the energy ratios with and without Coulomb interaction
noninteractingViral = kineticEnergy[0]/potentialEnergy[0]
interactingViral = kineticEnergy[1]/potentialEnergy[1]


# Plotting the energy ratios as function of the varying parameter omega
fig, ax = plt.subplots()
ax.plot(omega, noninteractingViral, color = '#BA8BCB', label='Non interacting')
ax.plot(omega, interactingViral, color = '#FEB144', label='Interacting')
ax.set_title(r'Energy ratio as function of $\omega$', fontsize=20)
ax.set_xlabel(r'$\omega$', fontsize=15)
ax.set_ylabel(r'$\frac{\langle T \rangle}{\langle V \rangle}$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/VirialTheorem.pdf')
