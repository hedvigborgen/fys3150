import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
style = ['--', '-']
colors = ['#1A7DA8', '#890C41']


# Defining arguments for C++ script
temperatures = [1, 2.4]
MCCs = 1000
try:
	L = sys.argv[1]
except:
	L = int(input('Enter size of matrix: '))


# Compiling the C++ script
#subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])
# Executing C++ script
os.system(f'./main.exe {L} 2 {MCCs}')


# Defining parameters
T_start = 2.0
T_end = 2.3
Ntemps = 100
h = (T_end - T_start)/(Ntemps - 1)


# Creating arrays for storing expectation values and temperatures
temp = np.linspace(T_start, T_end, 100)
expEnergy = np.zeros(len(temp))
expEnergySquared = np.zeros(len(temp))
expMagneticMoment = np.zeros(len(temp))
expMagneticMomentSquared = np.zeros(len(temp))


# Reading files written by C++ program
for i in range(Ntemps):
    T = T_start + i*h
    values = np.loadtxt(f'../output/{i}.dat')
    expEnergy[i] = values[-1, 1]/values[-1, 0]
    expEnergySquared[i] = values[-1, 2]/values[-1, 0]
    expMagneticMoment[i] = values[-1, 3]/values[-1, 0]
    expMagneticMomentSquared[i] = values[-1, 4]/values[-1, 0]


# Plotting the expectation value of the energy as function of temperature
fig, ax = plt.subplots()
ax.plot(temp, expEnergy, color = '#890C41')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Energy [J]', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/exp_energy_parallel_{L}_{MCCs}.pdf')


# Plotting the expectation value of the magnetic moment as function of temperature
fig, ax = plt.subplots()
ax.plot(temp, expMagneticMoment, color = '#890C41')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Magnetization', fontsize=15)
ax.set_title(f'Expectation value for the magnetic moment', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/exp_magnetization_parallel_{L}_{MCCs}.pdf')


# Numerical values of heat capacity and susceptibility
beta = 1/temp
cv = beta/temp*(expEnergySquared - expEnergy**2)
chi = beta*(expMagneticMomentSquared - expMagneticMoment**2)


# Plotting the numerical specific heat capacity as function of temperature
fig, ax = plt.subplots()
ax.plot(temp, cv, color='#1A7DA8', label=f'Numerical heat capacity')
ax.set_title(f'The specific heat capacity $C_V$ as a function of temperature $T$', fontsize=20)
ax.set_xlabel('$T$[J]', fontsize=15)
ax.set_ylabel('$C_V(T)$', fontsize=15)
fig.savefig(f'../plots/specific_heat_capacity_parallel_{L}_{MCCs}.pdf')


# Plotting the magnetic susceptibility as function of temperature
fig, ax = plt.subplots()
ax.plot(temp, chi, color='#1A7DA8', label=f'Numerical heat susceptibility')
ax.set_title(f'The susceptibility $\chi$ as a function of temperature $T$', fontsize=20)
ax.set_xlabel('$T$[J]', fontsize=15)
ax.set_ylabel('$\chi(T)$', fontsize=15)
fig.savefig(f'../plots/susceptibility_parallel_{L}_{MCCs}.pdf')
