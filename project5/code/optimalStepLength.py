import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

variations = 8
alpha0 = 0.3
deltaAlpha = 0.3
equilibrationTime = 100_000
MCCs = 1_000_000
charge = 1
whichMethod = 0
write = 'at the end'
steps = np.linspace(0.1, 3.5, 10)
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
            acceptedChanges[j-2,i] = float(vals[3])

    infile.close()


colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#F3E448', '#FC89B2', '#8ED3DB']

# Plotting the variation of the energy
# as function of MCCs for different values of alpha
fig, ax = plt.subplots()
for i, alpha_ in enumerate(alpha):
    ax.plot(steps, acceptedChanges[i]/MCCs * 100, color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)

ax.plot(steps, np.ones(len(steps))*50, ':', color = '#000000')
ax.set_title(r'Percentage of accepted changes as function of step length', fontsize=20)
ax.set_xlabel('Accepted changes', fontsize=15)
ax.set_ylabel(r'$\delta$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/optimalStepLength.pdf')