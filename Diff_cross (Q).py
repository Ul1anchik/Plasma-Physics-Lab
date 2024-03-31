import math
import matplotlib.pyplot as plt
import numpy as np


def poly_Legendre(filename):
    EjL = [0.5E+03, 0.5E+03, 0.15E+05, 0.5E+06, 0.15E+07, 0.3E+07, 0.35E+07, 0.4E+07, 0.7E+07]
    EjU = [0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09, 0.6357E+09]
    coef_cheb = np.loadtxt(filename)

    for j in range(coef_cheb.shape[0]):
        poly_cheb = np.polynomial.Chebyshev(coef_cheb[j])

        BOTTOM = EjL[j]
        TOP = EjU[j]
        E = BOTTOM
        N = 100
        STEP = pow(TOP / BOTTOM, 1 / N)

        coef_Legendre = []
        poly_Legendre = []

        while E <= TOP:
            Ej = 2 * (math.log(E / EjL[j]) / math.log(EjU[j] / EjL[j])) - 1
            coef_Leg = poly_cheb(Ej)
            coef_Legendre.append(coef_Leg)
            E *= STEP

        # while E > BOTTOM:
        x = np.linspace(-1, 1, 8)
        poly_leg = np.polynomial.legendre.legval(x, coef_Legendre, tensor=True)
        poly_Legendre.append(poly_leg)

    return poly_Legendre


def build_table(poly, filename): 
    EjL = 0.5000 * 10**3
    EjU = 0.6357 * 10**9
    
    BOTTOM = EjL
    TOP = EjU
    N = 100
    STEP = pow(TOP / BOTTOM, 1 / N)
    x = []
    y = []
    E = BOTTOM

    Legendre = poly_Legendre(filename)

    while E <= TOP:
        diff_cross_0 = diff_cross_sigma(poly, EjL, EjU, E)
        x.append(E)
        for i in range(len(Legendre)):
            diff_cross =  diff_cross_0* Legendre[i]
            y.append(diff_cross)
        E *= STEP
    return (x, y)

def diff_cross_sigma(poly, EjL, EjU, E):
    Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
    diff_cross = math.exp(poly(Ej))
    return diff_cross

def plot_diff_cross_section(ax: plt.Axes, filname, title, file):
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)

    x, y = build_table(poly, file)
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
        fig.add_subplot(1, 2, 1), "b_-1_p.txt", "Differential cross section $D(D,p) T$", 'coef_cheb_for_Leg_D(D,p)T.txt'
    )
    plot_diff_cross_section(
        fig.add_subplot(1, 2, 2), "b_-1_n.txt",  "Differential cross section $D(D,n)^3 He$", 'coef_cheb_for_Leg_D(D,n)T.txt'
    )
    plt.show()



  

