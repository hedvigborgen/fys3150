import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Defining input arguments for C++ & temperatures
part = 'c'
MCCs = 100_000
WhichMatrix = 2
L = 2
LL = L**2
temperature = np.linspace(0.5,5,1000)
temperature_num = np.linspace(0.5,5,10)

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'IsingModel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Arrays for storing expectation values
expEnergy = np.zeros(len(temperature_num))
expEnergySquared = np.zeros(len(temperature_num))
expMagneticMoment = np.zeros(len(temperature_num))
expMagneticMomentSquared = np.zeros(len(temperature_num))

# Looping over the different temperatures, running C++ script
for i, T in enumerate(temperature_num):
    os.system(f'./main.exe {part} {MCCs} {L} {T} {WhichMatrix}')
    values = np.loadtxt(f'../output/part4c_{T}.dat')
    MCCs_max = values[-1,0]
    expEnergy[i] = values[-1,1]
    expEnergySquared[i] = values[-1,2]
    expMagneticMoment[i] = values[-1,3]
    expMagneticMomentSquared[i] = values[-1,4]

    # For temperature T = 1 J
    if T == 1:
        # Calculating numerical values for heat capacity and susceptibility
        beta = 1/T
        cv_numerical = beta/T*(expEnergySquared[i] - expEnergy[i]**2)
        chi_numerical = beta*(expMagneticMomentSquared[i] - expMagneticMoment[i]**2)

        # Calculating analytical values
        e_exp = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))/LL
        e_squared_exp = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))/LL/LL
        m_exp = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))/LL
        m_squared_exp = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))/LL/LL
        cv_analytical = beta/T*(e_squared_exp - e_exp**2)
        chi_analytical = beta*(m_squared_exp - m_exp**2)

        # Printing numerical and analytical values for comparison
        print("For T = 1.0 J:")
        print("Numerical expectation value for the energy = ", expEnergy[i])
        print("Analytical expectation value for the energy = ", e_exp)
        print("Numerical expectation value for the energy squared = ", expEnergySquared[i])
        print("Analytical expectation value for the energy squared = ", e_squared_exp)
        print("Numerical expectation value for the magnetic moment = ", expMagneticMoment[i])
        print("Analytical expectation value for the magnetic moment = ", m_exp)
        print("Numerical expectation value for the magnetic moment squared = ", expMagneticMomentSquared[i])
        print("Analytical expectation value for the magnetic moment squared = ", m_squared_exp)
        print("Numerical expectation value for the specific heat capacity = ", cv_numerical)
        print("Analytical expectation value for the specific heat capacity = ", cv_analytical)
        print("Numerical expectation value for the susceptibility = ", chi_numerical)
        print("Analytical expectation value for the susceptibility = ", chi_analytical)

        # Storing expectation values per cycle for T = 1 J
        # & calculating numerical heat capacity and susceptibility
        # as function of MCCs
        MCCs_T1 = values[:,0]
        expEnergy_T1 = values[:,1]
        expEnergySquared_T1 = values[:,2]
        expMagneticMoment_T1 = values[:,3]
        expMagneticMomentSquared_T1 = values[:,4]
        cv_numerical_T1 = beta/T*(expEnergySquared_T1 - expEnergy_T1**2)
        chi_numerical_T1 = beta*(expMagneticMomentSquared_T1 - expMagneticMoment_T1**2)

        # Calculating corresponding analytical values
        e_exp_T1 = np.ones(len(MCCs_T1))*e_exp
        e_squared_exp_T1 = np.ones(len(MCCs_T1))*e_squared_exp
        m_exp_T1 = np.ones(len(MCCs_T1))*m_exp
        m_squared_exp_T1 = np.ones(len(MCCs_T1))*m_squared_exp
        cv_analytical_T1 = np.ones(len(MCCs_T1))*cv_analytical
        chi_analytical_T1 = np.ones(len(MCCs_T1))*chi_analytical

        # Plotting the numerical heat capacity for T = 1 J as function of MCCs
        # to compare with analytical value
        fig, ax = plt.subplots()
        ax.plot(MCCs_T1, cv_analytical_T1, color="#B1C084", label=f"Analytical heat capacity" )
        ax.plot(MCCs_T1, cv_numerical_T1, color='#1A7DA8', label=f"Numerical heat capacity")
        ax.set_title(f"The specific heat capacity $C_V$ as a function of MCCs", fontsize=20)
        ax.set_xlabel("MCCs", fontsize=15)
        ax.set_ylabel("$c_V(T)$", fontsize=15)
        ax.set_ylim(-0.01, 0.02)
        ax.legend(fontsize=15)
        fig.savefig(f'../plots/part_c/cv_functionofMCCs_cycle_{MCCs}.pdf')

        # Plotting the numerical susceptibility for T = 1 J as function of MCCs
        # to compare with analytical value
        fig, ax = plt.subplots()
        ax.plot(MCCs_T1, chi_analytical_T1, color="#B1C084", label=f"Analytical susceptibility" )
        ax.plot(MCCs_T1, chi_numerical_T1, color='#1A7DA8', label=f"Numerical susceptibility")
        ax.set_title(f"Susceptibility $\chi$ as a function of MCCs", fontsize=20)
        ax.set_xlabel("MCCs", fontsize=15)
        ax.set_ylabel("$\chi(T)$ [J$^{-1}$]", fontsize=15)
        ax.set_ylim(-0.01, 0.01)
        ax.legend(fontsize=15)
        fig.savefig(f'../plots/part_c/chi_functionofMCCs_cycle_{MCCs}.pdf')


# Calculating numerical values for heat capacity and susceptibility for random matrix
beta_num = 1/temperature_num
cv_numerical = beta_num/temperature_num*(expEnergySquared - expEnergy**2)
chi_numerical = beta_num*(expMagneticMomentSquared - expMagneticMoment**2)

# Calculating analytical values for the expectation values and partition function
beta = 1/temperature
e_exp = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))/LL
e_squared_exp = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))/LL/LL
m_exp = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))/LL
m_squared_exp = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))/LL/LL

# Calculating analytical values for heat capacity and susceptibility
cv_analytical = beta/temperature*(e_squared_exp - e_exp**2)
chi_analytical = beta*(m_squared_exp - m_exp**2)

# Plotting the numerical specific heat capacity as function of temperature
fig, ax = plt.subplots()
ax.plot(temperature, cv_analytical, color="#B1C084", label=f"Analytical heat capacity" )
ax.plot(temperature_num, cv_numerical, color='#1A7DA8', label=f"Numerical heat capacity")
ax.set_title(f"The specific heat capacity $C_V$ as a function of temperature $T$", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$C_V(T)$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_c/cv_functionofTemp_MCCs_{MCCs}.pdf')

# Plotting the magnetic susceptibility as function of temperature
fig, ax = plt.subplots()
ax.plot(temperature, chi_analytical, color="#B1C084", label=f"Analytical heat susceptibility" )
ax.plot(temperature_num, chi_numerical, color='#1A7DA8', label=f"Numerical heat susceptibility")
ax.set_title(f"The susceptibility $\chi$ as a function of temperature $T$", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$\chi(T)$ [J$^{-1}$]", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_c/chi_functionofTemp_MCCs_{MCCs}.pdf')

# Plotting expectation value for the magnetization as function of temperature
fig, ax = plt.subplots()
ax.plot(temperature, m_exp, color="#B1C084", label=f"Analytical expectation value")
ax.plot(temperature_num, expMagneticMoment, "--", color="#1A7DA8", label=f"Numerical expectation value")
ax.set_title(f"Analytical and numerical expectation value of magnetization", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$<M>$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_c/compare_expectation_magnetization{MCCs}.pdf')

# Plotting expectation value for energy squared as function of temperature
fig, ax = plt.subplots()
ax.plot(temperature, e_exp, color="#B1C084", label=f"Analytical expectation value")
ax.plot(temperature_num, expEnergy, "--", color="#1A7DA8", label=f"Numerical expectation value")
ax.set_title(f"Analytical and numerical expectation value of energy", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$<E>$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/part_c/compare_expectation_energy{MCCs}.pdf')
