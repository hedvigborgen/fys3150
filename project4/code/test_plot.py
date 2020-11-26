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
temperature = np.linspace(0.5,5,1000)
temperature_num = np.linspace(0.5,5,10)


expEnergy = np.zeros(len(temperature_num))
expEnergySquared = np.zeros(len(temperature_num))
expMagneticMoment = np.zeros(len(temperature_num))
expMagneticMomentSquared = np.zeros(len(temperature_num))

for i, T in enumerate(temperature_num):
    # Compiling the C++ script
    subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])
    os.system(f"./main.exe 20 2 {MCCs} {T}")
    values = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
    MCCs_max = values[:,0][-1]
    expEnergy[i] = values[:,1][-1]/MCCs_max
    expEnergySquared[i] = values[:,2][-1]/MCCs_max
    expMagneticMoment[i] = values[:,3][-1]/MCCs_max
    expMagneticMomentSquared[i] = values[:,4][-1]/MCCs_max

    if T == 1:
        beta = 1/T
        cv_numerical = beta/T*(expEnergySquared[i] - expEnergy[i]**2)
        chi_numerical = beta*(expMagneticMomentSquared[i] - expMagneticMoment[i]**2)

        e_exp = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))/LL
        e_squared_exp = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))/LL/LL
        m_exp = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))/LL
        m_squared_exp = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))/LL/LL
        cv_analytical = beta/T*(e_squared_exp - e_exp**2)
        chi_analytical = beta*(m_squared_exp - m_exp**2)

        print("expEnergy=", expEnergy[i])
        print("Analytical expEnergy=", e_exp)
        print("expEnergySquared=", expEnergySquared[i])
        print("Analytical expEnergySquared=", e_squared_exp)
        print("expMagneticMoment=", expMagneticMoment[i])
        print("Analytical expMagneticMoment=", m_exp)
        print("expMagneticMomentSquared=", expMagneticMomentSquared[i])
        print("Analytical expMagneticMomentSquared=", m_squared_exp)
        print("cv_numerical=", cv_numerical)
        print("Analytical cv=", cv_analytical)
        print("chi_numerical=", chi_numerical)
        print("Analytical chi=", chi_analytical)


        MCCs_T1 = values[:,0]
        expEnergy_T1 = values[:,1]/MCCs_T1
        expEnergySquared_T1 = values[:,2]/MCCs_T1
        expMagneticMoment_T1 = values[:,3]/MCCs_T1
        expMagneticMomentSquared_T1 = values[:,4]/MCCs_T1
        cv_numerical_T1 = beta/T*(expEnergySquared_T1 - expEnergy_T1**2)
        chi_numerical_T1 = beta*(expMagneticMomentSquared_T1 - expMagneticMoment_T1**2)

        e_exp_T1 = np.ones(len(MCCs_T1))*e_exp
        e_squared_exp_T1 = np.ones(len(MCCs_T1))*e_squared_exp
        m_exp_T1 = np.ones(len(MCCs_T1))*m_exp
        m_squared_exp_T1 = np.ones(len(MCCs_T1))*m_squared_exp
        cv_analytical_T1 = np.ones(len(MCCs_T1))*cv_analytical
        chi_analytical_T1 = np.ones(len(MCCs_T1))*chi_analytical


        fig, ax = plt.subplots()
        ax.plot(MCCs_T1, cv_analytical_T1, color="#B1C084", label=f"Analyticalheat capacity" )
        ax.plot(MCCs_T1, cv_numerical_T1, color='#1A7DA8', label=f"Numerical heat capacity")
        ax.set_title(f"The specific heat capacity $C_V$ as a function of MCCs", fontsize=20)
        ax.set_xlabel("MCCs", fontsize=15)
        ax.set_ylabel("$c_V(T)$", fontsize=15)
        ax.legend(fontsize=15)
        fig.savefig(f'../plots/specific_heat_capacity_cycle_{MCCs}.pdf')


beta = 1/temperature
beta_num = 1/temperature_num

# numerical values of heat capacity and susceptibility for random matrix
cv_numerical = beta_num/temperature_num*(expEnergySquared - expEnergy**2)
chi_numerical = beta_num*(expMagneticMomentSquared - expMagneticMoment**2)

# analytical values of the expectation values and partition function
e_exp = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))/LL
e_squared_exp = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))/LL/LL
m_exp = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))/LL
m_squared_exp = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))/LL/LL

# analytical values of heat capacity and susceptibiliy
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
fig.savefig(f'../plots/specific_heat_capacity{MCCs}.pdf')

# Plotting the magnetic susceptibility as function of temperature
fig, ax = plt.subplots()
ax.plot(temperature, chi_analytical, color="#B1C084", label=f"Analytical heat susceptibility" )
ax.plot(temperature_num, chi_numerical, color='#1A7DA8', label=f"Numerical heat susceptibility")
ax.set_title(f"The susceptibility $\chi$ as a function of temperature $T$", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$\chi(T)$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/susceptibility{MCCs}.pdf')

# Expectation magnetization
fig, ax = plt.subplots()
ax.plot(temperature, m_exp, color="#B1C084", label=f"Analytical")
ax.plot(temperature_num, expMagneticMoment, "--", color="#1A7DA8", label=f"Numerical")
ax.set_title(f"Analytical and numerical expectation value of magnetization", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$<M>$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/compare_expectation_magnetization{MCCs}.pdf')

# Expectation energy squared
fig, ax = plt.subplots()
ax.plot(temperature, e_exp, color="#B1C084", label=f"Analytical")
ax.plot(temperature_num, expEnergy, "--", color="#1A7DA8", label=f"Numerical")
ax.set_title(f"Analytical and numerical expectation value of energy", fontsize=20)
ax.set_xlabel("$T$[J]", fontsize=15)
ax.set_ylabel("$<E>$", fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/compare_expectation_energy{MCCs}.pdf')