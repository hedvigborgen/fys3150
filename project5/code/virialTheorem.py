import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

tasks = ['VirialwithoutInteraction', 'VirialwithInteraction']
MCCs = 1_000_000
omega0 = 0.01
deltaOmega = 0.01
size = 100

omega = np.linspace(omega0, deltaOmega*size, size)
kineticEnergy = np.zeros((len(tasks), size))
potentialEnergy = np.zeros((len(tasks), size))

# Compiling and executing the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
for i, task in enumerate(tasks):
    subprocess.call(['./main.exe', task, f'{MCCs}', f'{omega0}', f'{deltaOmega}'])
    filename = f'../output/Energy_{task}.dat'
    infile = open(filename, 'r')
    lines = infile.readlines()

    for j, line in enumerate(lines):
        vals = line.split()
        kineticEnergy[i,j] = float(vals[1])
        potentialEnergy[i,j] = float(vals[2])

    infile.close()


fig, ax = plt.subplots()

noninteractingViral = kineticEnergy[0]/potentialEnergy[0]
interactingViral = kineticEnergy[1]/potentialEnergy[1]

ax.plot(omega, noninteractingViral, color = '#BA8BCB', label='Non interacting')
ax.plot(omega, interactingViral, color = '#FEB144', label='Interacting')
ax.set_title(r'Energy ratio as function of $\omega$', fontsize=20)
ax.set_xlabel(r'$\omega$', fontsize=15)
ax.set_ylabel(r'$\frac{\langle T \rangle}{\langle V \rangle}$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/VirialTheorem.pdf')
