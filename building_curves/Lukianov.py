import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, exp 

def build_table(const):  
    x=[i for i in range(100, 100000, 100 )]
    y = []
    for i in x:
        values_y = 1000*(const/i)*(exp((-46)/math.sqrt(i)))
        y.append(values_y)
    return x,y

def plot_total_cross_section(ax: plt.Axes):

    x, y = build_table(300)
    ax.plot(x, y, color="green")
    # ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
   

fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1)
)

plt.show()
     