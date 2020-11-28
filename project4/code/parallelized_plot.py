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
part = 'f'
MCCs = 100_000
WhichMatrix = 2
Ls = [20, 40, 60, 80, 100]

# try:
# 	L = sys.argv[1]
# except:
# 	L = int(input('Enter size of matrix: '))

# Defining parameters
T_start = 2.0
T_end = 2.3
Ntemps = 100
h = (T_end - T_start)/(Ntemps - 1)

# Temperature array
temp = np.linspace(T_start, T_end, 100)

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'IsingModel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])


# Looping over different sizes of matrix
for L in Ls:
	# Executing C++ script
	os.system(f'./main.exe {part} {MCCs} {WhichMatrix} {L}')

	# Creating arrays for storing expectation values
	expEnergy = np.zeros(len(temp))
	expEnergySquared = np.zeros(len(temp))
	expMagneticMoment = np.zeros(len(temp))
	expMagneticMomentSquared = np.zeros(len(temp))

	# Reading files written by C++ program
	for i in range(len(temp)):
	    T = T_start + i*h
	    values = np.loadtxt(f'../output/part4f_{i}_L_{L}.dat')
	    expEnergy[i] = values[0]
	    expEnergySquared[i] = values[1]
	    expMagneticMoment[i] = values[2]
	    expMagneticMomentSquared[i] = values[3]

	# Plotting the expectation value of the energy as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expEnergy, color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Energy [J]', fontsize=15)
	ax.set_title(f'Expectation value for the energy, L = {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/exp_energy_parallel_{L}_{MCCs}.pdf')


	# Plotting the expectation value of the magnetic moment as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expMagneticMoment, color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Magnetization', fontsize=15)
	ax.set_title(f'Expectation value for the magnetic moment, L = {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/exp_magnetization_parallel_{L}_{MCCs}.pdf')


	# Numerical values of heat capacity and susceptibility
	beta = 1/temp
	cv = beta/temp*(expEnergySquared - expEnergy**2)
	chi = beta*(expMagneticMomentSquared - expMagneticMoment**2)


	# Plotting the numerical specific heat capacity as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, cv, color='#1A7DA8', label=f'Numerical heat capacity')
	ax.set_title(f'The specific heat capacity $C_V$ as a function of temperature $T$, L = {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$C_V(T)$', fontsize=15)
	fig.savefig(f'../plots/specific_heat_capacity_parallel_{L}_{MCCs}.pdf')


	# Plotting the magnetic susceptibility as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, chi, color='#1A7DA8', label=f'Numerical heat susceptibility')
	ax.set_title(f'The susceptibility $\chi$ as a function of temperature $T$, L = {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$\chi(T)$', fontsize=15)
	fig.savefig(f'../plots/susceptibility_parallel_{L}_{MCCs}.pdf')
