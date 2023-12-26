import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 
from scipy.integrate import quad
import scipy.integrate as spi

n = 10**(29)
y = []
x = [i for i in range(200, 100000, 100)]
for i in x:
    values_y = 7*10**(-10)*(n*n/i**(2/3))*exp((-4.25*10**3)/i**(1/3))
    y.append(values_y)

plt.plot(x, y) 
plt.grid(True) 
plt.yscale("log")
# plt.ylim([10**(-22), 10])
plt.xscale("log")
plt.show()  