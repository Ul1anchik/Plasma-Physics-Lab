import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 

def build_table(table_kin_energy, constants):
    with open(table_kin_energy, 'r') as file:
        array = file.readlines()
        x = [float(line.strip()) for line in array]

    y = []
    with open(constants, 'r') as file:
        Const = dict()
        for line in file:
            key, value = line.split()
            Const[key] = float(value)
    for i in x:       
        values_y = (Const['A2']/(1+(Const['A3']*i - Const['A4'])**2))/(i*(exp(Const['A1']/math.sqrt(i))-1))
        y.append(values_y)

    return x,y
    
def plot_total_cross_section(ax: plt.Axes, table_kin_energy, constants):

    x,y = build_table(table_kin_energy, constants)
    ax.plot(x, y, color="orange")
    # ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")

fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "kinetick_energy_D.txt", "constants D.txt"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "kinetick_energy_T.txt", "constants T.txt"
)
plt.show()