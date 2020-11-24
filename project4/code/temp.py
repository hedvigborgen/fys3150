import numpy as np
import matplotlib.pyplot as plt
import os
T_ = [1, 2.4]
I_ = [1,2]
L = 20
MCCs = int(1e5)

fig, ax = plt.subplots()
ax.set(xscale="log")
for T in T_:
	os.system(f"./main.exe 20 {T} 1 {MCCs}")
	order = np.loadtxt("../output/OrderedOrientation.dat")

	ax.plot(order[:,0], order[:,1], label=f"order {T}")

for T in T_:
	os.system(f"./main.exe 20 {T} 2 {MCCs}")
	random = np.loadtxt("../output/RandomOrientation.dat")
	ax.plot(random[:,0], random[:,1], label=f"random {T}")


ax.legend()
plt.show()
