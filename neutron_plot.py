import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["figure.dpi"] = 200
plt.rcParams["font.family"] = "Times New Roman"

k = 1
l = 5
j = 5.5

particle = "Neutron"

data1 = np.loadtxt(f'{particle}/{particle}_k={k},l={l},j={j:.1f}.txt', skiprows=1)

x1 = data1[:, 0]
y1 = data1[:, 1]
n = len(x1)
lst = []
for i in range(n):
    if y1[i] != 0.0:
        lst.append(i)
    else:
        continue

A = []
E = []
for jj in range(len(lst)):
    A.append(x1[lst[jj]])
    E.append(y1[lst[jj]])
if l == 0:
    L = 's'
elif l == 1:
    L = 'p'
elif l == 2:
    L = 'd'
elif l == 3:
    L = 'f'
elif l == 4:
    L = 'g'
elif l == 5:
    L = 'h'

if j == 0.5:
    J = '{1/2}'
elif j == 1.5:
    J = '{3/2}'
elif j == 2.5:
    J = '{5/2}'
elif j == 3.5:
    J = '{7/2}'
elif j == 4.5:
    J = '{9/2}'
elif j == 5.5:
    J = '{11/2}'


plt.figure(figsize=(6.5, 5.5))
plt.title(rf"Neutron, ${k}{L}_{J}$", fontsize = 20)
plt.plot(A, E, 'o', color = 'r',label = rf"${k}{L}_{J}$", markersize=3.0)
plt.ylim(-50, 10)
plt.xlabel('A (N+Z)')
plt.ylabel('E [MeV]')
plt.legend(fontsize = 20)

#print(j)
plt.grid()
plt.savefig(f"{particle}/{particle}_{k}{L}{j:.1f}.png")
# plt.show()