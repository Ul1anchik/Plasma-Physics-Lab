import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 
from scipy.integrate import quad

n = 10**(29)

def build_table_rate_coef(filname):
    
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)
    EjU = 0.6357 * 10**9
    EjL = 0.5000 * 10**3
    BOTTOM = EjL
    TOP = EjU
    N = 100
    STEP = pow(TOP / BOTTOM, 1 / N)
    cross_section = []
    E = BOTTOM
    while E <= TOP:
        Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
        total_cross = math.exp(poly(Ej))
        cross_section.append(total_cross)
        E *= STEP
    
    Vo = 2.18769126379E+08 #atomic unit of velocity
    md = 3.34358320*10**(-27) # deuteron mass(kg)
    M = md*md/(md+md) # in centre masses
    
    def integrand(y, T):
            return y*exp((-M*Vo*Vo/2*T)*y)
    
    x = [i for i in range(200, 1000099, 99)]
    y = []      
    for sigma in cross_section:
        for T in range(200,100000,1000):
            result, none = quad(integrand, 0, np.inf, args=(T) )
            R = (2/math.sqrt(np.pi)*Vo*((M*Vo*Vo/2*T)**(3/2)))*sigma*result
            y.append(R)
    return x,y

def build_table_Bosch(constants): #H.S.Bosch

    x = [i for i in range(200, 100000, 100)]
    y = []
    with open(constants, 'r') as file:
        Const = dict()
        for line in file:
            key, value = line.split()
            Const[key] = float(value)

    for i in x:
        Q = i / (1 - ( (i*(Const['C2']+i*(Const['C4']+i*Const['C6'])))/(1+i*(Const['C3']+i*(Const['C5']+i*Const['C7'])))))
        eps = ((Const['Bg']*Const['Bg'])/(4*Q))**(1/3)
        sigma = Const['C1'] * Q * math.sqrt(eps/(Const['mc']*Const['mc']*i**(3)))*exp(-3*eps)
        R = n*n*sigma
        y.append(R)

    return x,y

def build_table_Lukianov():
    y = []
    x = [i for i in range(200, 100000, 100)]
    for i in x:
        values_y = 7*10**(-10)*(n*n/i**(2/3))*exp((-4.25*10**3)/i**(1/3))
        y.append(values_y)
    return x,y

def plot_total_cross_section(ax: plt.Axes, filname, constants, title):

    x, y = build_table_rate_coef(filname)
    ax.plot(x, y, 'o', color="blue", markersize = 1, label='using a formula through a single integral')

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
    ax.set_xlim([10**(3), 10**9])
    ax.set_ylim([10**(-38), 10**16])


fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "bkj_Dp.txt", "building_rate_coef/constants_T_bosch.txt",  "Maxwellian rate coefficient $D(D,p) T$"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "bkj_Dn.txt", "building_rate_coef/constants_D_bosch.txt", "Maxwellian rate coefficient $D(D,n)^3 He$"
)
plt.show()
