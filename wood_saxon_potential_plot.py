import matplotlib.pyplot as plt
import numpy as np
plt.rcParams["figure.dpi"] = 200
plt.rcParams["font.family"] = "Times New Roman"

r0 = 1.27
a = 0.67
alpha = 1./137.036
c = 299792458 * 1.e15
hbar = 6.58211899 * 1.e-22
e2 = hbar * c * alpha
m_n = 938.272013 / c ** 2
effective_mass = m_n
pp = hbar ** 2 / (2 * effective_mass)

V_WS = lambda r, N, Z: (-51 + 33 * (N-Z)/(N+Z)) /(1 + np.exp((r-r0 * (N+Z) ** (1/3))/a))


def V_LS(r, N, Z, l, j, isospin):
    if isospin == 1:
        V0 = 0.44 * (51 + 33 * (N-Z)/(N+Z))
    elif isospin == -1:
        V0 = 0.44 * (51 - 33 * (N-Z)/(N+Z))
    else:
        print("Wrong isospin input!")
    R_q = r0 * (N+Z) ** (1./3.)
    t = (r - R_q)/a
    V_LS = -(r0 ** 2.) * np.exp(t) / a / ((1 + np.exp(t))) ** 2. / r * \
        (j * (j+1) - (l) * (l+1) - 3./4.) * V0 / 2.
    return V_LS


def V_C(r, N, Z, isospin):
    R_p = r0 * (N+Z) ** (1./3.)
    
    if isospin == -1:
        V_C = 0.
    
    elif isospin == 1:
        if r <= R_p:
            V_C = Z * e2 * (3 - (r/R_p) ** 2.) / (2. * R_p)
        
        elif r > R_p:
            V_C = Z * e2 / r
    else:
        print("Invalid Isospin!")
    return V_C    

def V_centrifugal(r, l):
    V_centrifugal = pp * ((l) * (l+1)) / (r**2.)
    return V_centrifugal



def V_effective(r, N, Z, l, j, isospin):
    V_effective = V_LS(r, N, Z, l, j, isospin) + V_C(r, N, Z, isospin) + \
                 V_WS(r, N, Z) + V_centrifugal(r, l)
    return V_effective





N = 10
Z = 10
k = 1
l = 1
j = 3/2
isospin = 1
x = np.arange(1.e-1, 10, 0.1)
n = len(x)
V_effective_array = []
for i in range(n):
    V = V_effective(x[i], N, Z, l, j, isospin)
    V_effective_array.append(V)
if isospin == 1:
    particle = 'Proton'
    color = 'b'
elif isospin == -1:
    particle = 'Neutron'
    color = 'r'

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
plt.plot(x, V_effective_array, color = color, label = f'N = {N}, Z = {Z} \n' +
                                        rf"${k}{L}_{J}$")
plt.title(f'Wood-Saxon Effective Potential of {particle}', fontsize=20)
plt.xlabel('r [fm]', fontsize=13)
plt.ylabel('V [MeV]', fontsize=13)
plt.ylim(-75., 50.)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc='upper right', fontsize=20)
plt.grid()
# plt.savefig(f'../../Desktop/research/meeting_0211/beamer/\
# Wood_Saxon_Potential{isospin}.png')
plt.savefig(f'../../Desktop/research/meeting_0211/beamer/\
Wood_Saxon_EffectivePotential{isospin}.png')

plt.show()