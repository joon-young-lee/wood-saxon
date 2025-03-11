import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl  
mpl.rc('font',family='Times New Roman')
plt.rcParams["figure.dpi"] = 300
plt.rcParams["font.family"] = "Times New Roman"

data = np.loadtxt('Radial_WF.txt', skiprows=1)
r_array = data[:, 0]
u = data[:, 1]
normalization = np.trapz(u, x=None, dx=r_array[3]-r_array[2])
u = u/normalization
plt.plot(r_array, u, color='b',label=r'$u = r\psi_\mathrm{Radial}$')
# plt.plot(r_array, u/r_array, color='r', label=r'$\psi_\mathrm{Radial}$')
plt.title(r"Proton $u$, $1p_{3/2}$", fontsize=20)
plt.xlabel('r [fm]', fontsize=15)
# plt.xticks(fontsize=18)
plt.legend(fontsize=20)
plt.grid()
plt.savefig('../../Desktop/research/meeting_0211/beamer/\
plot_u_1p1.5.png')
plt.show()




plt.plot(r_array, u/r_array, color='r', label=r'$\psi_\mathrm{Radial}$')
plt.title(r"Proton Radial Wave Function, $1p_{3/2}$", fontsize=20)
plt.xlabel('r [fm]', fontsize=15)

plt.legend(fontsize=20)
# plt.xticks(fontsize=18)
plt.grid()
plt.savefig('../../Desktop/research/meeting_0211/beamer/\
plot_psi_1p1.5.png')
plt.show()