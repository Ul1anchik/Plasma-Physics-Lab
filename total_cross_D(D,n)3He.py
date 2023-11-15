import math
import matplotlib.pyplot as plt
import numpy as np

EjL = 0.5000 * 10**3
EjU = 0.6357 * 10**9

def dzeta_i(E):
    return 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
    
coef = np.loadtxt("C:/Users/ryg20/OneDrive/Рабочий стол/физика плазмы/bkj_Dn.txt")
poly = np.polynomial.Chebyshev(coef)

BOTTOM = EjL
TOP = EjU

get_range = lambda: range(math.ceil(BOTTOM), math.floor(TOP), (math.floor(TOP) - math.ceil(BOTTOM)) // 100)


x = [i for i in get_range()]
y = []

for n in get_range():
    E = n
    Ej = dzeta_i(E)
    total_cross = math.exp(poly(Ej))
    y.append(total_cross)
    print(Ej, total_cross)

fig = plt.figure()
fig.suptitle('Total cross section D(D,n)3He')
ax = fig.add_subplot(1, 1, 1)
# ax.plot(x, y, color='black', alpha=1.4, marker='^', markevery=(0, 5), label='D(D,n)3He')
plt.plot(x, y, color = 'black')
plt.xlabel('$ \mathcal{E}$ (eV/u)')
plt.ylabel('$ \sigma$ (mb)')
plt.ylim([10**(-6), 10**2])


# plt.grid(True, color = "grey", linewidth = "1.0", linestyle = "-.")
ax.grid(which='minor', alpha=0.5, linestyle = "--" )
ax.grid(which='major', alpha=0.5, linestyle = "--")
# plt.loglog(x,y)
# plt.minorticks_on()
# plt.grid(True, which="both", axis="both", alpha=0.5, linestyle = "--")
# plt.show()

plt.xscale('log')
plt.yscale('log')

plt.show()
