import numpy as np
import matplotlib.pyplot as plt
import os

os.subprocess.run(['make compile'])

rho_max = np.linspace(3.0, 6.0, 31)
print(rho_max)
n = 100
for rho in rho_max:
    os.subprocess.run(['./main.exe', str(n), str(rho)])
