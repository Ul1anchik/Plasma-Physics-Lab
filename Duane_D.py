import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from scipy.optimize import curve_fit 
from numpy import array, exp 
import math
   
with open("kinetick_energy_D.txt", 'r') as file:
        array = file.readlines()
        values_x = [float(line.strip()) for line in array]

with open("cross_section D.txt", 'r') as file:
        array = file.readlines()
        values_y = [float(line.strip()) for line in array]

# with open("cross_section D_error.txt", 'r') as file:
#         array = file.readlines()
#         yerr = [float(line.strip()) for line in array]

A = 1527.0
B = 511000.0
C = 345.0E-09
D = 1.254

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
# plt.errorbar(values_x, values_y, xerr = 0, yerr=yerr, fmt='-')
plt.show() 