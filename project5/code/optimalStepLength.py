import numpy as np
import matplotlib.pyplot as plt
import subprocess
# from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import scipy

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

variations = 8
alpha0 = 0.3
deltaAlpha = 0.3
equilibrationTime = 100_000
MCCs = 1_000_000
charge = 1
whichMethod = 0
write = 'at the end'
steps = np.linspace(0.1, 3.5, 50)
alpha = np.linspace(0.2, 2.30, 8)


subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])


acceptedChanges = np.zeros((len(alpha), len(steps)))

# Compiling the C++ script
for i, step in enumerate(steps):
    subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{step}', f'{charge}', f'{whichMethod}', write])

    infile = open(f'../output/EnergyasFunctionofAlpha0_%.2f.dat' %step)
    lines = infile.readlines()

    for j, line in enumerate(lines):
        vals = line.split()
        if 2 <= j:
            acceptedChanges[j-2,i] = float(vals[3])/MCCs * 100

    infile.close()



distance = np.ones(len(alpha))*1000
idx = np.zeros(len(alpha), dtype='int')

for i in range(len(alpha)):
    for j in range(len(steps)):
        if distance[i] > np.abs(acceptedChanges[i,j] - 50):
            distance[i] = np.abs(acceptedChanges[i,j] - 50)
            idx[i] = j


# intercept1 = np.array([steps[idx[0]], acceptedChanges[0, idx[0]]])
# intercept2 = np.array([steps[idx[1]], acceptedChanges[1, idx[1]]])
# intercept3 = np.array([steps[idx[2]], acceptedChanges[2, idx[2]]])
# intercept4 = np.array([steps[idx[3]], acceptedChanges[3, idx[3]]])
# intercept5 = np.array([steps[idx[4]], acceptedChanges[4, idx[4]]])
# intercept6 = np.array([steps[idx[5]], acceptedChanges[5, idx[5]]])
# intercept7 = np.array([steps[idx[6]], acceptedChanges[6, idx[6]]])
# intercept8 = np.array([steps[idx[7]], acceptedChanges[7, idx[7]]])


# colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#FF6663','#F3E448', '#FC89B2', '#8ED3DB', '#1A7DA8']
#
# # Plotting the variation of the energy
# # as function of MCCs for different values of alpha
# fig, ax = plt.subplots()
# for i, alpha_ in enumerate(alpha):
#     ax.plot(steps, acceptedChanges[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
# ax.plot(steps, np.ones(len(steps))*50, ':', color = '#000000')
# ax.set_title(r'Percentage of accepted changes as function of step length', fontsize=20)
# ax.set_xlabel('Accepted changes', fontsize=15)
# ax.set_ylabel(r'$\delta$', fontsize=15)
# ax.legend(fontsize=15)
# fig.savefig(f'../plots/acceptedChanges.pdf')


# Finding best step length for
idealStep = np.zeros(len(alpha))
for i in range(len(alpha)):
    idealStep[i] = np.log(steps[idx[i]])

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(alpha, idealStep)
linearRegression = alpha*slope + intercept


fig, ax = plt.subplots()

<<<<<<< HEAD
ax.plot(alpha, linearRegression, color = '#BA8BCB', label=r'$\delta =$ %.3f$\alpha$ + %.3f' %(slope, intercept))
=======
ax.plot(alpha, linearRegression, color = '#BA8BCB')
>>>>>>> 63103c87d82d7beb8e815a585521eb301b49a333
for i in range(len(idealStep)):
    ax.plot(alpha[i], idealStep[i], 'o', color = '#000000')
ax.set_title(r'The optimal step length as function of alpha', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\ln{\delta}$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/optimalStepLength_log.pdf')
