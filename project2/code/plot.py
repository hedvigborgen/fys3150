import numpy as np
import matplotlib.pyplot as plt
import subprocess

subprocess.run(['make', 'compile'])
lines = 2
n = 100
rho_max = 5.0
omega_r = [0.01, 0.5, 1.0, 5.0]
for omega in omega_r:
    res = subprocess.run(['./main.exe', str(n), str(rho_max), str(omega)], capture_output=True, text=True)
    out = res.stdout.split(r'\n')
    print(out[:lines])
    with open(f'../output/res_{omega:.2f}.txt', 'w') as outfile:
        for line in out[lines:]:
            outfile.write(line)
