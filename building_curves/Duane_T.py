
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit 
from numpy import array, exp 
import math
   
with open("building_curves/kinetick_energy_T.txt", 'r') as file:
        array = file.readlines()
        values_x = [float(line.strip()) for line in array]

with open("building_curves/cross_section T.txt", 'r') as file:
        array = file.readlines()
        values_y = [float(line.strip()) for line in array]

A = 1457.7
B = 372000.0
C = 439.0E-09
D = 1.220

y_fit4 = []
for x in values_x:       
    y = (B/(1+(C*x - D)**2))/(x*(exp(A/math.sqrt(x))-1))
    y_fit4.append(y)

plt.plot(values_x, values_y, 'x', label="y - original") 
plt.plot(values_x, y_fit4) 
plt.xlabel('x') 
plt.ylabel('y') 
plt.legend(loc = 'best', fancybox = True, shadow = True) 
plt.grid(True) 
plt.yscale("log")
plt.xscale("log")
plt.show() 

