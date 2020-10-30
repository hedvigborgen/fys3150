import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

files = ["euler", "verlet"]

for i, e in enumerate(files):
    with open(f'../output/{e}.xyz', 'wb') as outfile:
        subprocess.call(['./main.exe', str(e)], stdout=outfile)

def f(filename, numTimesteps, numberOfBodies):
    infile = open(filename, 'r')
    lines = infile.readlines()
    x_Sun, y_Sun, z_Sun = np.zeros(n), np.zeros(n), np.zeros(n)
    x_Earth, y_Earth, z_Earth = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(0,numTimesteps,2):
        line = lines[i]
        values = line.split()
        x_Sun[i], y_Sun[i], z_Sun[i] = float(vals[0]), float(vals[0]), float(vals[0])
        print(x_Sun[i], y_Sun[i], z_Sun[i])
        input()
        line = lines[i+1]
        values = line.split()
        x_Earth[i], y_Earth[i], z_Earth[i] = float(vals[0]), float(vals[0]), float(vals[0])
        print(x_Earth[i], y_Earth[i], z_Earth[i])
        input()
    infile.close()
    return numerical, theoretical


