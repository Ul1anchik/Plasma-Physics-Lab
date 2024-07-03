import numpy as np

COEFS_FOR_LEGENDRE = "coef_cheb_for_Leg_D(D,n)He.txt"
COEFS_FOR_ZERO_ANGLE = "coef_cheb_for_zero_angle.txt"

E = 1.3e6


# IAEA
DIFF_CROSS_1_TABLE = {
    6e4: {"S0": 0.76, "A": [0.719, 0.279, 0.002]},
    0.15e6: {"S0": 3.39, "A": [0.633, 0.356, 0.011]},
    0.4e6: {"S0": 9.4, "A": [0.532, 0.426, 0.042]},
    0.65e6: {"S0": 14.3, "A": [0.474, 0.444, 0.082]},
    0.9e6: {"S0": 17.8, "A": [0.432, 0.443, 0.125]},
    1.3e6: {"S0": 21.9, "A": [0.381, 0.422, 0.189, 0.008]},
    1.8e6: {"S0": 25.2, "A": [0.334, 0.392, 0.251, 0.023]},
}


def diff_cross_1(s0, a, angle):
    return s0 * np.polynomial.legendre.legval(x=np.cos(angle), c=a)


# atomic and nuclear data
def non_linear_change(e, EjL=0.5e3, EjU=0.6357e9):
    return 2 * (np.log(e / EjL)) / (np.log(EjU / EjL)) - 1


ZERO_B = np.loadtxt(COEFS_FOR_ZERO_ANGLE)  # b_{k,-1} k = 0 .. 15
B = np.loadtxt(COEFS_FOR_LEGENDRE)  # b_{k,j} j = 0 .. 8 k = 0 .. 15


def s0_2(e):
    poly = np.polynomial.Chebyshev(ZERO_B)
    return np.exp(poly(non_linear_change(e)))


def diff_cross_2(e, angle):
    legendre_coef = []
    for j in range(8):
        legendre_coef.append(np.polynomial.Chebyshev(B[j])(non_linear_change(e)))
        legendre_coef.append(0)
    return s0_2(e) * np.polynomial.legendre.legval(x=np.cos(angle), c=legendre_coef)
