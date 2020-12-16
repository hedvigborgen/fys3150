import numpy as np
import matplotlib.pyplot as plt
import subprocess
from sklearn.linear_model import LinearRegression
#from scipy.optimize import curve_fit
import scipy

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#FF6663','#F3E448', '#FC89B2', '#8ED3DB', '#1A7DA8']

# Defining input parameters
task = 'StepLength'
MCCs = 1_000_000
whichMethod = 0
variations = 8
alpha0 = 0.2
deltaAlpha = 0.3
steps = np.linspace(0.1, 3.5, 50)
alpha = np.linspace(0.2, 2.30, 8)

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Array for storing number of accepted configurations as function of alpha
acceptedChanges = np.zeros((len(alpha), len(steps)))

# Executing the C++ script and reading in the data
for i, step in enumerate(steps):
    subprocess.call(['./main.exe', task, f'{MCCs}', f'{whichMethod}', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{step}'])

    infile = open(f'../output/EnergyasFunctionofAlpha0_%.2f.dat' %step)
    lines = infile.readlines()

    for j, line in enumerate(lines):
        vals = line.split()
        if 2 <= j:
            acceptedChanges[j-2,i] = float(vals[3])/MCCs * 100

    infile.close()


# Plotting the acceptance rate as function of step length
fig, ax = plt.subplots()
for i, alpha_ in enumerate(alpha):
    ax.plot(steps, acceptedChanges[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
ax.plot(steps, np.ones(len(steps))*50, ':', color = '#000000')
ax.set_title(r'Percentage of accepted changes as function of step length', fontsize=20)
ax.set_xlabel(r'$\delta$', fontsize=15)
ax.set_ylabel('Accepted changes / Total number of MCCs', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/acceptedChanges.pdf')


# Initializing arrays for finding the step length providing an acceptance rate of 50%
distance = np.ones(len(alpha))*1000
idx = np.zeros(len(alpha), dtype='int')

# Finding the step length providing an acceptance rate of 50%
for i in range(len(alpha)):
    for j in range(len(steps)):
        if distance[i] > np.abs(acceptedChanges[i,j] - 50):
            distance[i] = np.abs(acceptedChanges[i,j] - 50)
            idx[i] = j


# Array for storing the logaritm of the optimal step length for different alphas
idealStep = np.zeros(len(alpha))
for i in range(len(alpha)):
    idealStep[i] = np.log(steps[idx[i]])

# Trying to approximate the relation between the step length and alpha
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(alpha, idealStep)
linearRegression = alpha*slope + intercept

# Plotting the approximated relation between the step length and alpha
fig, ax = plt.subplots()
ax.plot(alpha, linearRegression, color = '#BA8BCB', label=r'$\delta = e^{%.3f\alpha + %.3f}$' %(slope, intercept))
for i in range(len(idealStep)):
    ax.plot(alpha[i], idealStep[i], 'o', color = '#000000')
ax.set_title(r'The optimal step length as function of alpha', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\ln{\delta}$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/optimalStepLength.pdf')
