
import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, exp 

def build_table(table_kin_energy, table_cross):

    with open(table_kin_energy, 'r') as file:
        array = file.readlines()
        x = [float(line.strip()) for line in array]
    with open(table_cross, 'r') as file:
        array = file.readlines()
        y = [float(line.strip()) for line in array]
    return x,y
    
def plot_total_cross_section(ax: plt.Axes, table_kin_energy, table_cross):

    x, y = build_table(table_kin_energy, table_cross)
    ax.plot(x, y, color="red")
    # ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")

fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "kinetick E bosch.txt", "bocsh T cross.txt"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "kinetick E bosch.txt", "bocsh D cross.txt"
)
plt.show()