import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["font.family"] = "Times New Roman"
x = np.arange(0, 15, 0.1)
r0 = 1.27
a = 0.67
f = lambda r, N, Z: (-51 + 33 * (N-Z)/(N+Z)) /(1 + np.exp((r-r0 * (N+Z) ** (1/3))/a))
N = 8
Z = 8
plt.plot(x, f(x, N, Z), label = f'N = {N}, Z = {Z}')
plt.title('Wood-Saxon Potential', fontsize=20)
plt.xlabel('r [fm]', fontsize=20)
plt.ylabel('V [MeV]', fontsize=20)
plt.legend(loc='lower right')
plt.savefig('Wood_Saxon_Potential.png')
plt.show()