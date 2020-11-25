import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

#for nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')



L = 2
LL = L**2
MCCs = 100000
temperature = 1

k_b = 1
beta = 1/temperature


# Compiling the C++ script
subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])
os.system(f"./main.exe 20 temperature 2 {MCCs}")
values = np.loadtxt("../output/RandomOrientation.dat")
cycle = values[:,0]
expEnergy = values[:,1]*LL
expEnergySquared = values[:,2]*LL*LL
expMagneticMoment = values[:,3]*LL
expMagneticMomentSquared = values[:,4]*LL*LL


#numerical values of heat capacity and susceptibiliy for random matrix
cv_numerical = beta/temperature*(expEnergySquared - expEnergy**2)
chi_numerical = beta*(expMagneticMomentSquared - expMagneticMoment**2)



#analytical values of the expectation values and partition function
e_squared_exp = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))
e_exp = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))
m_exp = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))
m_squared_exp = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))

#analytical values of heat capacity and susceptibiliy
cv_analytical = beta/temperature*(e_squared_exp - e_exp**2)
chi_analytical = beta*(m_squared_exp - m_exp**2)

cv_analyticalArray = np.ones(len(cycle))
cv_analyticalArray *= cv_analytical

chi_analyticalArray = np.ones(len(cycle))
chi_analyticalArray *= cv_analytical


fig, ax = plt.subplots()
ax.plot(cycle, cv_numerical, color="#B1C084", label=f"Numerical heat capacity" )
ax.plot(cycle, cv_analyticalArray, color='#1A7DA8', label=f"Analytical heat capacity")
ax.set_title(f"The specific heat capacity $C_V$ as a function of temperature $T$", fontsize=20)
ax.set_xlabel("MCCs", fontsize=15)
ax.set_ylabel("$C_V(T)$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/specific_heat_capacity.pdf')

fig, ax = plt.subplots()
ax.plot(cycle, chi_numerical, color="#B1C084", label=f"Numerical heat susceptibility" )
ax.plot(cycle, chi_analyticalArray, color='#1A7DA8', label=f"Analytical heat susceptibility")
ax.set_title(f"The susceptibility $\chi$ as a function of temperature $T$", fontsize=20)
ax.set_xlabel("MCCs", fontsize=15)
ax.set_ylabel("$\chi(T)$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/susceptibility.pdf')


# fig, ax = plt.subplots()
# # Need to add deltaE in order to calculate the numerical values for
# ax.plot(temperature, np.sqrt(m_exp_squared), color="#B1C084")
# ax.set_title(f"The expectation value of magnetization as a function of temperature $T$", fontsize=20)
# ax.set_xlabel("$T$[J]", fontsize=15)
# ax.set_ylabel("$<E>$", fontsize=15)
# fig.savefig(f'../plots/expectation_magnetization.pdf')
#
# fig, ax = plt.subplots()
# ax.plot(temperature, np.sqrt(e_exp_squared), color="#B1C084")
# ax.set_title(f"The expectation value of energy as a function of temperature $T$", fontsize=20)
# ax.set_xlabel("$T$[J]", fontsize=15)
# ax.set_ylabel("$<M>$", fontsize=15)
# fig.savefig(f'../plots/expectation_energy.pdf')
