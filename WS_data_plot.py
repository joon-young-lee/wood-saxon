import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["font.family"] = "Times New Roman"
data1 = np.loadtxt('potential.txt', skiprows=1)
l = 2
j = 3/2
T = 1
x1 = data1[:, 1]
y1 = data1[:, 0]
plt.title('A = 16')
plt.plot(x1, y1, label=f"l = {l}\n\
j = {j} \n\
T = {T}")
plt.xlim(0.0, 10.0)
ylim = 50.0
plt.ylim(-ylim, ylim)
plt.xlabel('r [fm]')
plt.ylabel('V [MeV]')
plt.legend(fontsize=20)
plt.savefig(f"l={l}_j={j}_T=1.png")
plt.show()