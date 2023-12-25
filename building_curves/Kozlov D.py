import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, exp 

e = 10**(-17)
def build(): #Building D cross section aprox formula
    
    x = [93.8227*i for i in range(100, 80000, 100)]
    y = []
    for i in x:
        values_y = 4.42286E+03*(0.481 + 0.0297*i)*(4.25/i)*( (exp(-3.2413/math.sqrt(i)))/(( i - 3.20 )**2 + 6.92) ) 
        y.append(values_y)
    return x, y

def build_table(table_kin_energy, table_cross): #building D cross experemental

    with open(table_kin_energy, 'r') as file:
        array = file.readlines()
        x = [93.8227*float(line.strip()) for line in array]
    with open(table_cross, 'r') as file:
        array = file.readlines()
        y = [4.42286*float(line.strip()) for line in array]
    return x,y

def build2(kinE):
    with open(kinE, 'r') as file:
        array = file.readlines()
        x = [93.8227*float(line.strip()) for line in array]
    y = []
    for i in x:
        values_y = 4.42286E+03*(0.481 + 0.0297*i)*(4.25/i)*( (exp(-3.2413/math.sqrt(i)))/(( i - 3.20 )**2 + 6.92) ) 
        y.append(values_y)
    return x, y


def plot_total_cross_section(ax: plt.Axes, table_kin_energy, table_cross, kinE):

    x, y = build_table(table_kin_energy, table_cross) #experimental data
    ax.plot(x, y, "x", color="red")

    x, y = build() #aprox formula with E in range
    ax.plot(x, y, color="blue")

    x, y = build2(kinE) #aprox formula with experemental E
    ax.plot(x, y, color="green")

    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")

fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1),"building_curves/kozlov_Ekin D.txt", "building_curves/kozlov_cross_section D.txt", "building_curves/kozlov_Ekin D.txt"
)
plt.show()
