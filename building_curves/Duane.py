import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 

def build_table(constants):

    x = [i for i in range(10000,5000000,1000) ]
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

def build(kinE, cross):
    with open(kinE, 'r') as file:
        array = file.readlines()
        values_x = [float(line.strip()) for line in array]

    with open(cross, 'r') as file:
        array = file.readlines()
        values_y = [float(line.strip()) for line in array]
    return values_x, values_y

    
def plot_total_cross_section(ax: plt.Axes,  constants, kinE, cross):

    x,y = build_table(constants)
    ax.plot(x, y, color="orange")
    x, y = build(kinE, cross)
    ax.plot(x, y,"x", color="blue")
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")

fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1),  "building_curves/constants D.txt" , "building_curves/kinetick_energy_D.txt", "building_curves/cross_section D.txt")
plot_total_cross_section(
    fig.add_subplot(1, 2, 2),  "building_curves/constants T.txt" , "building_curves/kinetick_energy_T.txt" , "building_curves/cross_section T.txt"
)
plt.show()