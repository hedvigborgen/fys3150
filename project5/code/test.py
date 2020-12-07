import numpy as np
import matplotlib.pyplot as plt

alpha = np.linspace(0.5, 1.5, 11)


E = 0.5*(1 - alpha**2) + 3*alpha
dE = -0.5*2*alpha + 3

# for i in range(len(alpha)):
#     if abs(dE[i] - E[i]) < abs(dE[i-1] - E[i-1]):
#         if abs(dE[i] - E[i]) < abs(dE[i+1] - E[i+1]):
plt.plot(alpha, E, label='E')
plt.plot(alpha, dE, label='dE')
plt.legend()
plt.show()
