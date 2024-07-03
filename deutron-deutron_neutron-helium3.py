import numpy as np

COEFS_FOR_LEGENDRE = "coef_cheb_for_Leg_D(D,n)He.txt"
COEFS_FOR_ZERO_ANGLE = "coef_cheb_for_zero_angle.txt"

E = 1.3e6


# IAEA
A = [0.381, 0.422, 0.189, 0.008]  # e = 1.3e6
S_0 = 21.9  # e = 1.3e6


def diff_cross_1(e, angle):
    global S_0, A
    return S_0 * np.polynomial.legendre.legval(x=np.cos(angle), c=A)


# atomic and nuclear data
def non_linear_change(e):
    EjL = 0.5e3
    EjU = 0.6357e9
    return 2 * (np.log(e / EjL)) / (np.log(EjU / EjL)) - 1


ZERO_B = np.loadtxt(COEFS_FOR_ZERO_ANGLE)  # b_{k,-1} k = 0 .. 15
B = np.loadtxt(COEFS_FOR_LEGENDRE)  # b_{k,j} j = 0 .. 8 k = 0 .. 15


def s0_2(e):
    return np.exp(np.polynomial.Chebyshev(ZERO_B)(non_linear_change(e)))


def diff_cross_2(e, angle):
    legendre_coef = [
        np.polynomial.Chebyshev(B[j])(non_linear_change(e)) for j in range(8)
    ]
    return s0_2(e) * np.polynomial.legendre.legval(x=np.cos(angle), c=legendre_coef)


if __name__ == "__main__":
    for angle in np.linspace(0, np.pi, 10):
        print(diff_cross_1(E, angle), diff_cross_2(E, angle))
