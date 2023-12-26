import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp 
from scipy.integrate import quad
import scipy.integrate as spi


# coef = np.loadtxt(filname)
# poly = np.polynomial.Chebyshev(coef)
# EjU = 0.6357 * 10**9
# EjL = 0.5000 * 10**3
# BOTTOM = EjL
# TOP = EjU
# N = 100
# STEP = pow(TOP / BOTTOM, 1 / N)
# cross_section = []
# E = BOTTOM
# while E <= TOP:
#     Ej = 2 * (math.log(E / EjL) / math.log(EjU / EjL)) - 1
#     total_cross = math.exp(poly(Ej))
#     cross_section.append(total_cross)
#     E *= STEP
  

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

def build_table_rate_coef(filname):
    
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
    
    Vo = 2.18769126379E+08 #atomic unit of velocity
    md = 3.34358320*10**(-27) # deuteron mass(kg)
    M = md*md/(md+md) # in centre masses

    def integrand(y, T):
            return y*exp((-M*Vo*Vo/2*T)*y)
    
    x = [i for i in range(200, 1000099, 99)]
    y = []      
    for sigma in cross_section:
        for T in range(200,100000,1000):
            result, none = quad(integrand, 0, np.inf, args=(T) )
            R = (2/math.sqrt(np.pi)*((M*Vo*Vo/2*T)**(3/2)))*Vo*sigma*result
            y.append(R)
    return x,y


x, y = build_table_rate_coef("bkj_Dp.txt")
plt.scatter(x, y)
# plt.plot(x, y, 'o', color="blue", markersize = 1)
plt.xlabel("T")
plt.ylabel("R")
plt.xscale("log")
plt.yscale("log")
plt.show()

    
    
    