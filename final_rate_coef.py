import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp
from math import pi
from scipy.integrate import quad


n = 10**(20)
def build_table_rate_coef(filname):
    EjL = 0.5 * 10**3  # eV
    EjU = 0.6357 * 10**9  # eV
    V0 = 2.18769126379 * 10**6  # atomic unit of velocity (m/s)
    md = 3.34358320 / 1.661  # deuteron mass (a.e.m)
    mu = md**2 / (2 * md)  # in centre masses (a.e.m)
    
    A = mu * (V0**2) / 2000  # eV (a.e.m * m^2 / s^2)
    
    def gamma(T):
        return A / T
    
    def Ej(y):
        return 2 * math.log(A * y / EjL, EjU / EjL) - 1
    
    def integrand(y, T):
        poly_res = poly(Ej(y))
        sigma = exp(poly_res)
        f = y * exp(-y * gamma(T))
        return sigma * f
    
    coef = np.loadtxt(filname)
    poly = np.polynomial.Polynomial(coef)
    
    LOWER = 10**3
    UPPER = 10**9
    STEP = 1.5
    
    T = LOWER
    
    Ts = []
    Rs = []
    while T < UPPER:
        CONST = 2 / math.sqrt(pi) * V0 * gamma(T) ** (3 / 2)
        integral_res, *_ = quad(lambda y: integrand(y, T), 0, np.inf)
        R = CONST * integral_res
        Ts.append(T)
        Rs.append(integral_res)
        T *= STEP
    return Ts, Rs

def build_table_Bosch(constants): #H.S.Bosch

    x = [i for i in range(1, 100, 1)] #keV
    y = []
    with open(constants, 'r') as file:
        Const = dict()
        for line in file:
            key, value = line.split()
            Const[key] = float(value)
    
    for i in x:
        Q = i / (1 - ( (i*(Const['C2']+i*(Const['C4']+i*Const['C6'])))/(1+i*(Const['C3']+i*(Const['C5']+i*Const['C7'])))))
        eps = ((Const['Bg']*Const['Bg'])/(4*Q))**(1/3)
        sigma = Const['C1'] * Q * math.sqrt(eps/(Const['mc']*i**(3)))*exp(-3*eps)
        R = n*n*1/2*sigma
        y.append(sigma)

    return x,y


def build_table_Lukianov():
    y = []
    x = [i for i in range(11_000_000, 100_000_000, 1000)] #K
    for i in x:
        values_y = 14*10**(-10)*(1/i**(2/3))*exp((-4.25*10**3)/i**(1/3))
        y.append(values_y)
    return [8.617 * 10**(-8) * i for i in x],y


def plot_total_cross_section(ax: plt.Axes, filname, constants, title):

    x, y = build_table_rate_coef(filname)
    ax.plot(x, y, color="blue", markersize = 1, label='using a formula through a single integral')

    x, y = build_table_Bosch(constants) #by H.S.Bosch
    ax.plot(x, y, color="red", label='using H.S.Bosch formula') 

    x,y = build_table_Lukianov() #by Lukianov
    ax.plot(x, y, color="green", label='using S.U.Lukianov formula') 
    

    ax.set_title(title)
    ax.set_xlabel("T")
    ax.set_ylabel("R")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_ylim([10**(-38), 10**19])


fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "bkj_Dp.txt","building_rate_coef/constants_T_bosch2.txt",  "Maxwellian rate coefficient $D(D,p) T$"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "bkj_Dn.txt",  "building_rate_coef/constants_D_bosch2.txt", "Maxwellian rate coefficient $D(D,n)^3 He$"
)
plt.show()


