import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 
from scipy.integrate import quad
import scipy.integrate as spi

def build_table1(poly, filname):

    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)

    EjU = 0.6357 * 10**9
    EjL = 0.5000 * 10**3
    BOTTOM = EjL
    TOP = EjU
    N = 100
    STEP = pow(TOP / BOTTOM, 1 / N)
    cross_section = []
    E = BOTTOM
    while E <= TOP:
        Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
        total_cross = math.exp(poly(Ej))
        cross_section.append(total_cross)
        E *= STEP
    return cross_section

# T = [i for i in range(200, 100000, 100)]
Vo = 2.18769126379E+08 #atomic unit of velocity
md = 3.34358320*10**(-27) # deuteron mass(kg)
M = md*md/(md+md) # in centre masses
# gamma = M*Vo*Vo/2*T #dimensionless parameter
# integrative = (2/math.sqrt(np.pi)*Vo*((M*Vo*Vo/2*T)**(3/2)))*sigma*y*exp((-M*Vo*Vo/2*T)*y)


# # print(cross_section)
# sigma = 3
# T = 200
# result, none = quad(integrand, 0, np.inf - 0.001 , args=(T) )
# R = (2/math.sqrt(np.pi)*Vo*((M*Vo*Vo/2*T)**(3/2)))*sigma*result


# def integrand(y, sigma, T):
#         return (2/math.sqrt(np.pi)*Vo*((M*Vo*Vo/2*T)**(3/2)))*sigma*y*exp((-M*Vo*Vo/2*T)*y)
# y = []
# for sigma in cross_section:
#     for T in range(1000,100000,1000):
#         result, none = quad(integrand, 0, np.inf - 0.001 , args=(sigma, T) )
#         R = result
#         y.append(R)
# x = [i for i in range(1000, 1000900, 100)]


def integrand(y, T):
        return y*exp((-M*Vo*Vo/2*T)*y)
y = []      
for sigma in cross_section:
    for T in range(200,100000,1000):
        result, none = quad(integrand, 0, np.inf - 0.001 , args=(T) )
        R = (2/math.sqrt(np.pi)*Vo*((M*Vo*Vo/2*T)**(3/2)))*sigma*result
        y.append(R)

x = [i for i in range(200, 1000099, 99)]
plt.plot(x, y) 
plt.grid(True) 
plt.yscale("log")
plt.xscale("log")
plt.show() 

def plot_total_cross_section(ax: plt.Axes, filname, title):
    coef = np.loadtxt(filname)
    poly = np.polynomial.Chebyshev(coef)
    
    
    