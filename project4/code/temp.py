import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

temperature = [1, 2.4]
whichMatrix = [1,2]
L = 20
MCCs = int(1e5)

# try:
# 	L, MCCs = sys.argv[1:]
# except:
# 	L = int(input('Enter size of matrix: '))
# 	MCCs = int(input('Enter number of MMCs: '))


# Compiling the C++ script
#subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])

style = ["--", "-"]
colors = ['#1A7DA8', '#890C41']

# Plotting the mean energy as function of the number of MCCs,
#for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
ax.set(xscale="log")
for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 1 {MCCs}")
	ordered = np.loadtxt(f"../output/OrderedOrientation_{T}.dat")
	ax.plot(ordered[:,0], ordered[:,1], f"{style[0]}", color = colors[i], label=f"Ordered spins, T = {T}")


for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	ax.plot(random[:,0], random[:,1], f"{style[1]}", color = colors[i], label=f"Random spins, T = {T}")


ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Energy', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueEnergy.pdf')

# Plotting the magnetization (absolute values) as function of the number of MCCs,
#for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
ax.set(xscale="log")
for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 1 {MCCs}")
	ordered = np.loadtxt(f"../output/OrderedOrientation_{T}.dat")
	ax.plot(ordered[:,0], ordered[:,3], f"{style[0]}", color = colors[i], label=f"Ordered spins, T = {T}")


for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	ax.plot(random[:,0], random[:,3], f"{style[1]}", color = colors[i], label=f"Random spins, T = {T}")


ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Magnetization', fontsize=15)
ax.set_title(f'Expectation value for the magnetization (absolute value)', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueMagnetization.pdf')

# Plotting the total number of accepted configurations as function of the number of MCCs,
#for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
ax.set(xscale="log")
for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 1 {MCCs}")
	ordered = np.loadtxt(f"../output/OrderedOrientation_{T}.dat")
	ax.plot(ordered[:,0], ordered[:,5], f"{style[0]}", color = colors[i], label=f"Ordered spins, T = {T}")


for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	ax.plot(random[:,0], random[:,5], f"{style[1]}", color = colors[i], label=f"Random spins, T = {T}")


ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Number of flips', fontsize=15)
ax.set_title(f'Total number of accepted configurations (flips)', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/numberOfFlips.pdf')
