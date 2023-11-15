import math
import matplotlib.pyplot as plt
import numpy as np

EjL = 0.5000 * 10**3
EjU = 0.6357 * 10**9


def dzeta_i(E):
    return 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1


coef = np.loadtxt("bkj_Dp.txt")
poly = np.polynomial.Chebyshev(coef)

BOTTOM = EjL
TOP = EjU

get_range = lambda: range(
    math.ceil(BOTTOM), math.floor(TOP), (math.floor(TOP) - math.ceil(BOTTOM)) // 1000000
)


x = [i for i in get_range()]
y = []

for n in get_range():
    E = n
    Ej = dzeta_i(E)
    total_cross = math.exp(poly(Ej))
    y.append(total_cross)
    print(Ej, total_cross)

fig = plt.figure()
fig.suptitle("Total cross  D(D,p)3He")
ax = fig.add_subplot(1, 1, 1)

plt.plot(x, y, color="black", label="D(D,n)3He")
plt.xlabel("$ \mathcal{E}$ (eV/u)")
plt.ylabel("$ \sigma$ (mb)")

ax.grid(which="minor", alpha=0.5, linestyle="--")
ax.grid(which="major", alpha=0.5, linestyle="--")

plt.xscale("log")
plt.yscale("log")

plt.show()
