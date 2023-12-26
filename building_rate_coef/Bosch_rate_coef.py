import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 
from scipy.integrate import quad
import scipy.integrate as spi

n = 10**(29)

def build_table(constants): #H.S.Bosch

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

def plot_total_cross_section(ax: plt.Axes, constants, title):

    x,y = build_table(constants) 
    ax.plot(x, y, color="orange") 
    
    ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    


fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "building_rate_coef/constants_T_bosch.txt",  "Maxwellian rate coefficient $D(D,p) T$"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "building_rate_coef/constants_D_bosch.txt", "Maxwellian rate coefficient $D(D,n)^3 He$"
)
plt.show()
