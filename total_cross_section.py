import math
import matplotlib.pyplot as plt
import numpy as np


EjL = 0.5000 * 10**3
EjU = 0.6357 * 10**9

BOTTOM = EjL
TOP = EjU
N = 100
STEP = pow(TOP / BOTTOM, 1 / N)


def build_table(poly):
    x = []
    y = []
    E = BOTTOM

    while E <= TOP:
        Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
        total_cross = math.exp(poly(Ej))
        x.append(E)
        y.append(total_cross)
        print(Ej, total_cross)

        E *= STEP

    return (x, y)


def plot_total_cross_section(ax: plt.Axes, filname, title):
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)

    x, y = build_table(poly)
    ax.plot(x, y, color="black")
    ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")


fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "bkj_Dp.txt", "Total cross section $D(D,p) T$"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "bkj_Dn.txt", "Total cross section $D(D,n)^3 He$"
)
plt.show()
