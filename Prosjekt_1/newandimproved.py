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
h = (x[-1]-x[0])/(n-1)

a = -1
b = 2
c = -1


f = h**2*f(x)

# Forward substitution:
for i in range(1, n):
    f[i] += f[i-1]*i/(i + 1) # Oppdaterer til f_tilde

f[-1] /= (n+1)/n

# Backward substitution:
for j in range(n-2, -1, -1):
    f[j] = (f[j] + f[j+1])*(j)/(j+1) # Oppdaterer til u


#f[-1] /=

plt.plot(x, exact(x), label='Exact U(x)')
plt.plot(x, f, label='Numerical approximation')
plt.legend()
plt.xlabel('x')
plt.ylabel('U(x)')
plt.show()
