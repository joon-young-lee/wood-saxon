import matplotlib.pyplot as plt
import numpy as np
data1 = np.loadtxt('potential.txt', skiprows=1)

x1 = data1[:, 1]
y1 = data1[:, 0]

plt.plot(x1, -y1, label="n = 1")

plt.legend()
plt.savefig("Particle_in_a_box.png")
plt.show()