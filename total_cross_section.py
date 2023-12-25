import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 

def build_table1(poly): 
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
        Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
        total_cross = math.exp(poly(Ej))
        x.append(E)
        y.append(total_cross)
        # print(Ej, total_cross)

        E *= STEP

    return (x, y)

def build_table2(table_kin_energy, table_cross): #H.S.Bosch

    with open(table_kin_energy, 'r') as file:
        array = file.readlines()
        x = [float(line.strip()) for line in array]
    with open(table_cross, 'r') as file:
        array = file.readlines()
        y = [float(line.strip()) for line in array]
    return x,y

def build_table3(constants): #B.H.Duane

    x = [i for i in range(1000,5000000,1000)]
    y = []
    with open(constants, 'r') as file:
        Const = dict()
        for line in file:
            key, value = line.split()
            Const[key] = float(value)
    for i in x:       
        values_y = (Const['A2']/(1+(Const['A3']*i - Const['A4'])**2))/(i*(exp(Const['A1']/math.sqrt(i))-1))*1000
        y.append(values_y)
    return x,y

def build_table4(kin_energy, cross): #by Kozlov

    with open(kin_energy, 'r') as file:
        array = file.readlines()
        x = [93.8227*float(line.strip()) for line in array]
    with open(cross, 'r') as file:
        array = file.readlines()
        y = [4.42286*float(line.strip()) for line in array]
    return x,y

def build_table5(const = 300):  #by Lukianov
    x=[i for i in range(100, 100000, 100 )]
    y = []
    for i in x:
        values_y = 1000*(const/i)*(exp((-46)/math.sqrt(i)))
        y.append(values_y)
    return x,y

def plot_total_cross_section(ax: plt.Axes, filname, table_kin_energy, table_cross, constants, kin_energy, cross, title):
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)

    x, y = build_table1(poly)
    ax.plot(x, y, color="black")

    x, y = build_table2(table_kin_energy, table_cross) #by H.S.Bosch
    ax.plot(x, y, color="red") 

    x,y = build_table3(constants) #by B.H.Duane
    ax.plot(x, y, color="orange") 

    x,y = build_table4(kin_energy, cross) #by Kozlov
    ax.plot(x, y, color="blue") 

    x,y = build_table5(const = 300) #by Lukianov
    ax.plot(x, y, color="green") 
    
    ax.set_title(title)
    ax.set_xlabel("$ \mathcal{E}$ (eV/u)")
    ax.set_ylabel("$ \sigma$ (mb)")
    ax.grid(which="minor", alpha=0.5, linestyle="--")
    ax.grid(which="major", alpha=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([10**(-4), 1000])


fig = plt.figure()
plot_total_cross_section(
    fig.add_subplot(1, 2, 1), "bkj_Dp.txt", "building_curves/kinetick E bosch.txt", "building_curves/bocsh T cross.txt",  "building_curves/constants T.txt", "building_curves/kozlov_Ekin T.txt", "building_curves/kozlov_cross_section T.txt",  "Total cross section $D(D,p) T$"
)
plot_total_cross_section(
    fig.add_subplot(1, 2, 2), "bkj_Dn.txt", "building_curves/kinetick E bosch.txt", "building_curves/bocsh D cross.txt", "building_curves/constants D.txt", "building_curves/kozlov_Ekin D.txt", "building_curves/kozlov_cross_section D.txt", "Total cross section $D(D,n)^3 He$"
)
plt.show()
