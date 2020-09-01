import numpy as np
import matplotlib.pyplot as plt

# U(x):
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

# f(x):
def f(x):
    return 100*np.exp(-10*x)


n = 100
x = np.linspace(0,1,n)
h = (x[-1]-x[0])/n

a = -1
b = 2
c = -1


f = h**2*f(x)

# Forward substitution:
for i in range(n-1):
    f[i+1] += f[i]*(i + 1)/(i + 2) # Oppdaterer til f_tilde

f[-1] /= (n+1)/n

# Backward substitution:
for j in range(n-1, 0, -1):
    f[j-1] += f[j]*(j+1)/(j+2) # Oppdaterer til u


plt.plot(x, exact(x), label='Exact U(x)')
plt.plot(x, f, label='Numerical approximation')
plt.legend()
plt.xlabel('x')
plt.ylabel('U(x)')
plt.show()
