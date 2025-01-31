import matplotlib.pyplot as plt
import numpy as np
data1 = np.loadtxt('output1.txt', skiprows=1)
data2 = np.loadtxt('output2.txt', skiprows=1)
data3 = np.loadtxt('output3.txt', skiprows=1)
data4 = np.loadtxt('output4.txt', skiprows=1)
x1 = data1[:, 0]
y1 = data1[:, 1]
x2 = data2[:, 0]
y2 = data2[:, 1]
x3 = data3[:, 0]
y3 = data3[:, 1]
x4 = data4[:, 0]
y4 = data4[:, 1]
plt.plot(x1, y1, label="n = 1")
plt.plot(x2, y2, label="n = 2")
plt.plot(x3, y3, label="n = 3")
plt.plot(x4, y4, label="n = 4")
plt.legend()
plt.savefig("Particle_in_a_box.png")
plt.show()