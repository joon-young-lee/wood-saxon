import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["figure.dpi"] = 200
plt.rcParams["font.family"] = "Times New Roman"

P_data0 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=0,j=0.5.txt', skiprows=1)
P_data1 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=1,j=0.5.txt', skiprows=1)
P_data2 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=1,j=1.5.txt', skiprows=1)
P_data4 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=2,j=1.5.txt', skiprows=1)
P_data5 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=2,j=2.5.txt', skiprows=1)
P_data6 = np.loadtxt('./wood-saxon/Proton/Proton_k=1,l=3,j=2.5.txt', skiprows=1)

P_A0 = P_data0[:, 0]
P_E0 = P_data0[:, 1]
P_A1 = P_data1[:, 0]
P_E1 = P_data1[:, 1]
P_A2 = P_data2[:, 0]
P_E2 = P_data2[:, 1]
P_A4 = P_data4[:, 0]
P_E4 = P_data4[:, 1]
P_A5 = P_data5[:, 0]
P_E5 = P_data5[:, 1]
P_A6 = P_data6[:, 0]
P_E6 = P_data6[:, 1]
print(P_A0)
plt.scatter(P_A0, P_E0, label=r'$1s_{1/2}$')
plt.scatter(P_A1, P_E1, label=r'$1p_{1/2}$')
plt.scatter(P_A2, P_E2, label=r'$1p_{3/2}$')
plt.scatter(P_A4, P_E4, label=r'$1d_{3/2}$')
plt.scatter(P_A5, P_E5, label=r'$1d_{5/2}$')
plt.scatter(P_A6, P_E6, label=r'$1f_{5/2}$')
plt.title('Proton Energy of Orbitals', fontsize=20)
plt.xlabel(r'A (N+Z)', fontsize=18)
plt.ylabel(r'$E_n$ [MeV]', fontsize=18)
plt.ylim(-35., -2.)
plt.legend()
plt.savefig('../../Desktop/research/meeting_0211/beamer/P_orbital.png')
plt.show()