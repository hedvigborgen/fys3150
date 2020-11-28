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

# Defining parameters
part = "c"
MCCs = 100_000
WhichMatrix = 2
L = 2
LL = L*L

T_numerical = np.linspace(0.5, 5, 10)
T_analytical = np.linspace(0.5, 5, 1000)


# Creating arrays for storing numerical expectation values
expEnergy_numerical = np.zeros((len(T_numerical), MCCs))
expEnergySquared_numerical = np.zeros((len(T_numerical), MCCs))
expMagneticMoment_numerical = np.zeros((len(T_numerical), MCCs))
expMagneticMomentSquared_numerical = np.zeros((len(T_numerical), MCCs))

# Creating arrays for storing analytic expectation values
expEnergy_analytical = np.zeros((len(T_numerical), MCCs))
expEnergySquared_analytical = np.zeros((len(T_numerical), MCCs))
expMagneticMoment_analytical = np.zeros((len(T_numerical), MCCs))
expMagneticMomentSquared_analytical = np.zeros((len(T_numerical), MCCs))



# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Reading in expectation values for a random spin matrix with method 1
for i, T in enumerate(T_numerical):
    os.system(f'./main.exe {part} {MCCs} {WhichMatrix} {L} {T}')
    values = np.loadtxt(f'../output/part4c_{T}.dat')
    expEnergy_numerical[i,:] = values[:,1]
    expEnergySquared_numerical[i,:] = values[:,2]
    expMagneticMoment_numerical[i,:] = values[:,3]
    expMagneticMomentSquared_numerical[i,:] = values[:,4]

    if T == 1:
        # Calculating heat capacity and susceptibility for numerical results
        beta = 1/T
        cv_numerical = beta/T*(expEnergySquared_numerical[i] - expEnergy_numerical[i]**2)
        chi_numerical = beta*(expMagneticMomentSquared_numerical[i] - expMagneticMoment_numerical[i]**2)

        expEnergy_analytical = -8*np.sinh(8*beta)/(3+np.cosh(8*beta))/LL
        expEnergySquared_analytical = 8**2*(np.cosh(8*beta)/(3+np.cosh(8*beta)))/LL/LL
        expMagneticMoment_analytical = 2*(2+np.exp(8*beta))/(3+np.cosh(8*beta))/LL
        expMagneticMomentSquared_analytical = 8*(1+np.exp(8*beta))/(3+np.cosh(8*beta))/LL/LL
        cv_analytical = beta/T*(expEnergySquared_analytical - expEnergy_analytical**2)
        chi_analytical = beta*(expMagneticMomentSquared_analytical - expMagneticMoment_analytical**2)

        print("expEnergy = ", expEnergy_numerical[i])
        print("Analytical expEnergy = ", expEnergy_analytical)
        print("expEnergySquared = ", expEnergySquared_numerical[i])
        print("Analytical expEnergySquared = ", expEnergySquared_analytical)
        print("expMagneticMoment = ", expMagneticMoment_numerical[i])
        print("Analytical expMagneticMoment = ", expMagneticMoment_analytical)
        print("expMagneticMomentSquared = ", expMagneticMomentSquared_numerical[i])
        print("Analytical expMagneticMomentSquared = ", expMagneticMomentSquared_analytical)
        print("cv_numerical = ", cv_numerical)
        print("Analytical cv = ", cv_analytical)
        print("chi_numerical = ", chi_numerical)
        print("Analytical chi = ", chi_analytical)


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
        ax.plot(MCCs_T1[1:], cv_analytical_T1[1:], color="#B1C084", label=f"Analytical heat capacity" )
        ax.plot(MCCs_T1[1:], cv_numerical_T1[1:], color='#1A7DA8', label=f"Numerical heat capacity")
        ax.set_title(f"The specific heat capacity $C_V$ as a function of MCCs", fontsize=20)
        ax.set_xlabel("MCCs", fontsize=15)
        ax.set_ylabel("$c_V(T)$", fontsize=15)
        ax.legend(fontsize=15)
        fig.savefig(f'../plots/cv_functionofMCCs_cycle_{MCCs}.pdf')

        fig, ax = plt.subplots()
        ax.plot(MCCs_T1, chi_analytical_T1, color="#B1C084", label=f"Analytical susceptibility" )
        ax.plot(MCCs_T1, chi_numerical_T1, color='#1A7DA8', label=f"Numerical susceptibility")
        ax.set_title(f"Susceptibility $\chi$ as a function of MCCs", fontsize=20)
        ax.set_xlabel("MCCs", fontsize=15)
        ax.set_ylabel("$\chi(T)$", fontsize=15)
        ax.legend(fontsize=15)
        fig.savefig(f'../plots/chi_functionofMCCs_cycle_{MCCs}.pdf')
