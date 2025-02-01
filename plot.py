import matplotlib.pyplot as plt
import numpy as np
data1 = np.loadtxt('k=0,l=0,j=0.5.txt', skiprows=1)

x1 = data1[:, 0]
y1 = data1[:, 1]

plt.plot(x1, y1, 'o')
plt.ylim(-40, 10)
plt.legend()
plt.savefig("state.png")
plt.show()