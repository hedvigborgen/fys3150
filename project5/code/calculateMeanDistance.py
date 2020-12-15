import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']


# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

def meanDistance(whichMethod, step, alpha, beta, omega):
    task = 'Parameters'
    MCCs = 1_000_000

    meanDistance = np.zeros(len(omega))
    for i, omega_ in enumerate(omega):
        subprocess.call(['./main.exe', task, f'{MCCs}', f'{whichMethod}', f'{alpha}', f'{beta}', f'{omega_}'])

        infile = open(f'../output/EnergyasFunctionofParameters_{whichMethod}_%.2f_%.2f_%.2f.dat' %(alpha, beta, omega_))
        lines = infile.readlines()
        vals = lines[2].split()
        meanDistance[i] = float(vals[5])
        infile.close()
    return meanDistance


whichMethod = [0, 1]
alpha = np.array([1.0, 0.85])
step = np.exp(-0.518*alpha + 0.982)
beta = 1.0
omega = [0.01, 0.5, 1.0]

meanDistance0 = meanDistance(whichMethod[0], step[0], alpha[0], beta, omega) # Without interaction
meanDistance1 = meanDistance(whichMethod[1], step[1], alpha[1], beta, omega) # With interaction

print('Without Coulomb interaction:')
print(r'For omega = 0.001, <r_{1,2}> = %.2f.' %meanDistance0[0])
print(r'For omega = 0.5, <r_{1,2}> = %.2f.' %meanDistance0[1])
print(r'For omega = 1.0, <r_{1,2}> = %.2f.' %meanDistance0[2])

print('With Coulomb interaction:')
print(r'For omega = 0.001, <r_{1,2}> = %.2f.' %meanDistance1[0])
print(r'For omega = 0.5, <r_{1,2}> = %.2f.' %meanDistance1[1])
print(r'For omega = 1.0, <r_{1,2}> = %.2f.' %meanDistance1[2])
