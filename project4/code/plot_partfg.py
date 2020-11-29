import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Defining arguments for C++ script
part = 'f'
MCCs = 1_000_000
WhichMatrix = 2
Ls = [20, 40, 60, 80, 100]

# Defining parameters
T_start = 2.0
T_end = 2.3
Ntemps = 100
h = (T_end - T_start)/(Ntemps - 1)

# Temperature array
temp = np.linspace(T_start, T_end, 100)

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'IsingModel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Creating arrays for storing expectation values
expEnergy = np.zeros((len(Ls), len(temp)))
expEnergySquared = np.zeros((len(Ls), len(temp)))
expMagneticMoment = np.zeros((len(Ls), len(temp)))
expMagneticMomentSquared = np.zeros((len(Ls), len(temp)))

# Looping over different sizes of matrix
for i, L in enumerate(Ls):
	# Executing C++ script
	os.system(f'./main.exe {part} {MCCs} {L} {WhichMatrix}')

	# Reading output files from C++ program
	for j in range(len(temp)):
	    T = T_start + i*h
	    values = np.loadtxt(f'../output/part4f_{i}_L_{L}.dat')
	    expEnergy[i, j] = values[0]
	    expEnergySquared[i, j] = values[1]
	    expMagneticMoment[i, j] = values[2]
	    expMagneticMomentSquared[i, j] = values[3]

	# Plotting the expectation value of the energy as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expEnergy, color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Energy [J]', fontsize=15)
	ax.set_title(f'Expectation value for the energy, L = {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/part_fg/exp_energy_parallel_{L}_{MCCs}.pdf')

	# Plotting the expectation value of the magnetic moment as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expMagneticMoment, color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Magnetization', fontsize=15)
	ax.set_title(f'Expectation value for the magnetic moment, L = {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/part_fg/exp_magnetization_parallel_{L}_{MCCs}.pdf')

	# Calculating the numerical heat capacity and susceptibility
	beta = 1/temp
	cv = beta/temp*(expEnergySquared - expEnergy**2)
	chi = beta*(expMagneticMomentSquared - expMagneticMoment**2)

	# Plotting the numerical specific heat capacity as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, cv, color='#1A7DA8')
	ax.set_title(f'The specific heat capacity $C_V$ as a function of temperature $T$, L = {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$C_V(T)$', fontsize=15)
	fig.savefig(f'../plots/part_fg/specific_heat_capacity_parallel_{L}_{MCCs}.pdf')

	# Plotting the magnetic susceptibility as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, chi, color='#1A7DA8')
	ax.set_title(f'The susceptibility $\chi$ as a function of temperature $T$, L = {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$\chi(T)$', fontsize=15)
	fig.savefig(f'../plots/part_fg/susceptibility_parallel_{L}_{MCCs}.pdf')


# Estimating the curie temperature T_C:
# Calculating the numerical heat capacity and susceptibility for different Ls
beta = 1/temp

# L = 40
cv_40 = beta/temp*(expEnergySquared[1] - expEnergy[1]**2)
chi_40 = beta*(expMagneticMomentSquared[1] - expMagneticMoment[1]**2)
# L = 60
cv_60 = beta/temp*(expEnergySquared[2] - expEnergy[2]**2)
chi_60 = beta*(expMagneticMomentSquared[2] - expMagneticMoment[2]**2)
# L = 80
cv_80 = beta/temp*(expEnergySquared[3] - expEnergy[3]**2)
chi_80 = beta*(expMagneticMomentSquared[3] - expMagneticMoment[3]**2)
# L = 100
cv_100 = beta/temp*(expEnergySquared[4] - expEnergy[4]**2)
chi_100 = beta*(expMagneticMomentSquared[4] - expMagneticMoment[4]**2)


# Estimating T_c from the heat capacity
# L = 40
cv_max_40 = np.max(cv_40)
index_40 = np.where(cv_max_40 == np.max(cv_40))[0]
T_c_40_cv = temp[index_40]
# L = 60
cv_max_60 = np.max(cv_60)
index_60 = np.where(cv_max_60 == np.max(cv_46))[0]
T_c_60_cv = temp[index_60]
# L = 80
cv_max_80 = np.max(cv_80)
index_80 = np.where(cv_max_80 == np.max(cv_80))[0]
T_c_80_cv = temp[index_80]
# L = 100
cv_max_100 = np.max(cv_100)
index_100 = np.where(cv_max_100 == np.max(cv_100))[0]
T_c_100_cv = temp[index_100]


# Estimating T_c from the susceptibility
# L = 40
chi_max_40 = np.max(chi_40)
index_40 = np.where(chi_max_40 == np.max(chi_40))[0]
T_c_40_chi = temp[index_40]
# L = 60
chi_max_60 = np.max(chi_60)
index_60 = np.where(chi_max_60 == np.max(chi_46))[0]
T_c_60_chi = temp[index_60]
# L = 80
chi_max_80 = np.max(chi_80)
index_80 = np.where(chi_max_80 == np.max(chi_80))[0]
T_c_80_chi = temp[index_80]
# L = 100
chi_max_100 = np.max(chi_100)
index_100 = np.where(chi_max_100 == np.max(chi_100))[0]
T_c_100_chi = temp[index_100]


# Printing the resulting values for T_c
print('Estimated T_c from specific heat capacity:')
print('For L = 40, T_c = %.5f' %T_c_40_cv)
print('For L = 60, T_c = %.5f' %T_c_60_cv)
print('For L = 80, T_c = %.5f' %T_c_80_cv)
print('For L = 100, T_c = %.5f' %T_c_100_cv)
print('Estimated T_c from magnetic susceptibility:')
print('For L = 40, T_c = %.5f' %T_c_40_chi)
print('For L = 60, T_c = %.5f' %T_c_60_chi)
print('For L = 80, T_c = %.5f' %T_c_80_chi)
print('For L = 100, T_c = %.5f' %T_c_100_chi)


# Plotting the numerical specific heat capacity as function of temperature for different Ls
fig, ax = plt.subplots()
ax.plot(temp, cv_40, 'o', color='#FEB144', label=f'L = 40')
ax.plot(temp, cv_60, 'o', color='#9EE09E', label=f'L = 60')
ax.plot(temp, cv_80, 'o', color='#1A7DA8', label=f'L = 80')
ax.plot(temp, cv_100, 'o', color='#FF6663', label=f'L = 100')
ax.set_title(f'The specific heat capacity $C_V$ for different lattice sizes $L$', fontsize=20)
ax.set_xlabel('$T$ [J]', fontsize=15)
ax.set_ylabel('$C_V(T)$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_fg/estimatingTc_cv.pdf')


# Plotting the magnetic susceptibility as function of temperature for different Ls
fig, ax = plt.subplots()
ax.plot(temp, chi_40, 'o', color='#FEB144', label=f'L = 40')
ax.plot(temp, chi_60, 'o', color='#9EE09E', label=f'L = 60')
ax.plot(temp, chi_80, 'o', color='#1A7DA8', label=f'L = 80')
ax.plot(temp, chi_100, 'o', color='#FF6663', label=f'L = 100')
ax.set_title(f'The susceptibility $\chi$ for different lattice sizes L', fontsize=20)
ax.set_xlabel('$T$ [J]', fontsize=15)
ax.set_ylabel('$\chi(T)$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_fg/estimatingTc_chi.pdf')
