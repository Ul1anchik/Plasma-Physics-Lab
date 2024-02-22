import math
import matplotlib.pyplot as plt
import numpy as np

def build_table(poly): 
    EjL = 0.5000 * 10**3
    EjU = 0.6357 * 10**9
    
    BOTTOM = EjL
    TOP = EjU
    N = 100
    STEP = pow(TOP / BOTTOM, 1 / N)
    x = []
    y = []
    E = BOTTOM

    while E <= TOP:
        diff_cross = diff_cross_sigma(poly, EjL, EjU, E)
        x.append(E)
        y.append(diff_cross)
        # print(Ej, total_cross)
        E *= STEP
    return (x, y)

def diff_cross_sigma(poly, EjL, EjU, E):
    Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
    diff_cross = math.exp(poly(Ej))
    return diff_cross

def plot_diff_cross_section(ax: plt.Axes, filname, title):
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)

    x, y = build_table(poly)
    ax.plot(x, y, color="black")
    
    ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ d\sigma/d\Omega$ (mb/sr)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")

if __name__ == "__main__":
    fig = plt.figure()
    plot_diff_cross_section(
        fig.add_subplot(1, 2, 1), "b_-1_p.txt", "Differential cross section $D(D,p) T$"
    )
    plot_diff_cross_section(
        fig.add_subplot(1, 2, 2), "b_-1_n.txt",  "Differential cross section $D(D,n)^3 He$"
    )
    plt.show()
