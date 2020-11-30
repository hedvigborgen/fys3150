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
T_start = 2.15
T_end = 2.45
Ntemps = 20
h = (T_end - T_start)/(Ntemps - 1)

# Temperature array
temp = np.linspace(T_start, T_end, Ntemps)

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
	    T = T_start + j*h
	    values = np.loadtxt(f'../output/part4f_{j}_L_{L}.dat')
	    expEnergy[i, j] = values[0]
	    expEnergySquared[i, j] = values[1]
	    expMagneticMoment[i, j] = values[2]
	    expMagneticMomentSquared[i, j] = values[3]

	# Plotting the expectation value of the energy as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expEnergy[i], color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Energy [J]', fontsize=15)
	ax.set_title(f'Expectation value for the energy, $L =$ {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/part_fg/exp_energy_parallel_{L}_{MCCs}.pdf')

	# Plotting the expectation value of the magnetic moment as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, expMagneticMoment[i], color = '#890C41')
	ax.tick_params(axis='both', which='major', labelsize=15)
	ax.set_xlabel(r'$T$ [J]', fontsize=15)
	ax.set_ylabel(r'Magnetization', fontsize=15)
	ax.set_title(f'Expectation value for the magnetic moment, $L =$ {L}', fontsize=20)
	fig.tight_layout()
	fig.savefig(f'../plots/part_fg/exp_magnetization_parallel_{L}_{MCCs}.pdf')

	# Calculating the numerical heat capacity and susceptibility
	beta = 1/temp
	cv = beta**2*(expEnergySquared[i] - expEnergy[i]**2)
	chi = beta*(expMagneticMomentSquared[i] - expMagneticMoment[i]**2)

	# Plotting the numerical specific heat capacity as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, cv, color='#1A7DA8')
	ax.set_title(f'The specific heat capacity $C_V$ as a function of temperature $T$, $L =$ {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$C_V(T)$', fontsize=15)
	fig.savefig(f'../plots/part_fg/specific_heat_capacity_parallel_{L}_{MCCs}.pdf')

	# Plotting the magnetic susceptibility as function of temperature
	fig, ax = plt.subplots()
	ax.plot(temp, chi, color='#1A7DA8')
	ax.set_title(f'The susceptibility $\chi$ as a function of temperature $T$, $L =$ {L}', fontsize=20)
	ax.set_xlabel('$T$ [J]', fontsize=15)
	ax.set_ylabel('$\chi(T)$ [J$^{-1}$]', fontsize=15)
	fig.savefig(f'../plots/part_fg/susceptibility_parallel_{L}_{MCCs}.pdf')


# Estimating the curie temperature T_C:
# Calculating the numerical heat capacity and susceptibility for different Ls
beta = 1/temp

# L = 20
cv_20 = beta/temp*(expEnergySquared[0] - expEnergy[0]**2)
chi_20 = beta*(expMagneticMomentSquared[0] - expMagneticMoment[0]**2)
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
# L = 20
index_20 = np.where(cv_20 == np.max(cv_20))[0][0]
T_c_20_cv = temp[index_20]
# L = 40
index_40 = np.where(cv_40 == np.max(cv_40))[0][0]
T_c_40_cv = temp[index_40]
# L = 60
index_60 = np.where(cv_60 == np.max(cv_60))[0][0]
T_c_60_cv = temp[index_60]
# L = 80
index_80 = np.where(cv_80 == np.max(cv_80))[0][0]
T_c_80_cv = temp[index_80]
# L = 100
index_100 = np.where(cv_100 == np.max(cv_100))[0][0]
T_c_100_cv = temp[index_100]


# Estimating T_c from the susceptibility
# L = 20
index_20 = np.where(chi_20 == np.max(chi_20))[0][0]
T_c_20_chi = temp[index_20]
# L = 40
index_40 = np.where(chi_40 == np.max(chi_40))[0][0]
T_c_40_chi = temp[index_40]
# L = 60
index_60 = np.where(chi_60 == np.max(chi_60))[0][0]
T_c_60_chi = temp[index_60]
# L = 80
index_80 = np.where(chi_80 == np.max(chi_80))[0][0]
T_c_80_chi = temp[index_80]
# L = 100
index_100 = np.where(chi_100 == np.max(chi_100))[0][0]
T_c_100_chi = temp[index_100]


# Printing the resulting values for T_c(L)
print('Estimated T_c from specific heat capacity:')
print('T_c(20) = %.5f' %T_c_20_cv)
print('T_c(40) = %.5f' %T_c_40_cv)
print('T_c(60) = %.5f' %T_c_60_cv)
print('T_c(80) = %.5f' %T_c_80_cv)
print('T_c(100) = %.5f' %T_c_100_cv)
print('Estimated T_c from magnetic susceptibility:')
print('T_c(20) = %.5f' %T_c_20_cv)
print('T_c(40) = %.5f' %T_c_40_chi)
print('T_c(60) = %.5f' %T_c_60_chi)
print('T_c(80) = %.5f' %T_c_80_chi)
print('T_c(100) = %.5f' %T_c_100_chi)


# Plotting the numerical specific heat capacity as function of temperature for different Ls
fig, ax = plt.subplots()
ax.plot(temp, cv_20, 'o', color='#BA8BCB', label=f'L = 20')
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
ax.plot(temp, chi_20, 'o', color='#BA8BCB', label=f'L = 20')
ax.plot(temp, chi_40, 'o', color='#FEB144', label=f'L = 40')
ax.plot(temp, chi_60, 'o', color='#9EE09E', label=f'L = 60')
ax.plot(temp, chi_80, 'o', color='#1A7DA8', label=f'L = 80')
ax.plot(temp, chi_100, 'o', color='#FF6663', label=f'L = 100')
ax.set_title(f'The susceptibility $\chi$ for different lattice sizes $L$', fontsize=20)
ax.set_xlabel('$T$ [J]', fontsize=15)
ax.set_ylabel('$\chi(T)$ [J$^{-1}$]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_fg/estimatingTc_chi.pdf')


# Estimating T_c(L = infinity)
from sklearn.linear_model import LinearRegression
import scipy

T_c_40 = T_c_40_chi
T_c_60 = T_c_60_chi
T_c_80 = T_c_80_chi
T_c_100 = T_c_100_chi

T_c = np.array([T_c_40, T_c_60, T_c_80, T_c_100])
L = np.array(Ls[1:])
x = 1/L

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, T_c)
print('Estimated T_c(L = infinity) = %.5f, with standard deviation = %.5f' %(intercept, std_err))

# Plotting the linear regression of T_c
line = slope*x + intercept
fig, ax = plt.subplots()
ax.plot(x[0], line[0], 'o', color='#FEB144', label=f'$T_c(L=\infty) = %.5f$' %intercept)
ax.plot(x, line, color='#FF6663')
ax.set_title(f'Linear regression of the critical temperature', fontsize=20)
ax.set_xlabel('$1/L$', fontsize=15)
ax.set_ylabel('$T_C$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_fg/estimatingTc.pdf')
