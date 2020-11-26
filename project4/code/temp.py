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
MCCs = int(10000)

# try:
# 	L, MCCs = sys.argv[1:]
# except:
# 	L = int(input('Enter size of matrix: '))
# 	MCCs = int(input('Enter number of MMCs: '))


# Compiling the C++ script
subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])

style = ["--", "-"]
colors = ['#1A7DA8', '#890C41']

# # Plotting the mean energy as function of the number of MCCs,
# # for both ordered and random matrices with T = 1 and 2.4.
#
fig, ax = plt.subplots()
ax.set(xscale="log")
for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 1 {MCCs}")
	ordered = np.loadtxt(f"../output/OrderedOrientation_{T}.dat")
	MCCs_ordered = ordered[:,0]
	expEnergy_ordered = ordered[:,1]#/(MCCs_ordered)
	# print(MCCs_ordered)
	# print(expEnergy_ordered)
    # expEnergySquared[i] = values[:,2][-1]/MCCs_max
    # expMagneticMoment[i] = values[:,3][-1]/MCCs_max
    # expMagneticMomentSquared[i] = values[:,4][-1]/MCCs_max
	ax.plot(MCCs_ordered, expEnergy_ordered, f"{style[0]}", color = colors[i], label=f"Ordered spins, T = {T}")


for i, T in enumerate(temperature):
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	MCCs_random = random[:,0]
	expEnergy_random = random[:,1]#/(MCCs_random)
	# print(MCCs_random)
	# print(expEnergy_random)
	ax.plot(MCCs_random, expEnergy_random, f"{style[1]}", color = colors[i], label=f"Random spins, T = {T}")


ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Energy', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueEnergy{MCCs}.pdf')

#Plotting the magnetization (absolute values) as function of the number of MCCs,
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
fig.savefig(f'../plots/expecationvalueMagnetization{MCCs}.pdf')

# Plotting the total number of accepted configurations as function of the number of MCCs,
# for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
for i, T in enumerate(temperature):
	ax.set(xscale="log")
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	ax.plot(random[:,0], random[:,5], f"{style[1]}", color = colors[i], label=f"T = {T}")

ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Number of flips', fontsize=15)
ax.set_title(f'Total number of accepted flips', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/numberOfFlipsT_{T}_{MCCs}.pdf')


# Plotting accepted flips as function of temperature in percent
temperature_array = np.linspace(0,5,100)
numberOfFlips = np.zeros(len(temperature_array))

for i, T in enumerate(temperature_array):
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	values = np.loadtxt(f"../output/RandomOrientation_{T}.dat")
	flips = values[:,5]
	numberOfFlips[i] = flips[-1]

fig, ax = plt.subplots()
ax.plot(temperature_array, numberOfFlips/max(numberOfFlips)*100, "o",color = '#890C41')
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'$T$[J]', fontsize=15)
ax.set_ylabel('Accepted flips [$\%$]', fontsize=15)
ax.set_title(f'Accepted flips as function of temperature', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/flipsOfTemperature{MCCs}.pdf')



#Indices where temperature is between 2 and 2.5
index1 = np.where(temperature_array>2.5)[0][0]
index2 = np.where(temperature_array<2)[0][-1]

#Calulating slope of number of flips
grad = np.gradient(numberOfFlips)
#Max derivative = max slope
maxgrad = np.amax(abs(grad[index2:index1]))

#Index of temp with highest slope rate
index_t = np.where(abs(grad) == maxgrad)[0][0]
#Numerical curie temperature
t_c = temperature_array[index_t]

fig, ax = plt.subplots()
ax.plot(temperature_array, numberOfFlips/max(numberOfFlips)*100, "o", color = '#890C41')
ax.plot(t_c, numberOfFlips[index_t]/max(numberOfFlips)*100, "o", color = "#4F6DD5", label = '$T_C=$ %.2f'%t_c)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'$T$[J]', fontsize=15)
ax.set_ylabel('Accepted flips [$\%$]', fontsize=15)
ax.set_title(f'Accepted flips as function of temperature', fontsize=20)
ax.legend(fontsize=15)
fig.tight_layout()
fig.savefig(f'../plots/flipsOfTemperature{MCCs}.pdf')


#P(E)
for i, T in enumerate(temperature):
	fig, ax = plt.subplots()
	os.system(f'./main.exe 20 {T} 2 {MCCs}')
	values = np.loadtxt(f'../output/RandomOrientation_{T}.dat')
	expEnergy = values[:,1][-1]
	expEnergySquared = values[:,2][-1]
	totalEnergy = values[:,6]
	variance = np.abs(expEnergySquared - expEnergy**2)
	print(variance)
	ax.hist(totalEnergy, bins='auto', color= '#890C41')
	ax.set_xlabel(r'$E$[J]', fontsize=15)
	ax.set_ylabel('$P(E)$', fontsize=15)
	ax.set_title('Probability of each energy state, $P(E)$', fontsize=20)
	fig.savefig(f'../plots/probability_{T}_{MCCs}.pdf')