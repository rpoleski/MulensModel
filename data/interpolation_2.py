import math
import numpy as np
from math import sin, cos, sqrt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.special import ellipk, ellipe
# These are complete elliptic integrals of the first and the second kind.
from sympy.functions.special.elliptic_integrals import elliptic_pi as ellip3

def _get_magnification_WM94(u, rho):
        if u == rho:
            u2 = u**2
            a = np.pi / 2. + np.arcsin((u2 - 1.) / (u2 + 1.))
            return (2./u + (1.+u2) * a / u2) / np.pi

        a_1 = 0.5 * (u + rho) * (4. + (u-rho)**2)**.5 / rho**2
        a_2 = -(u - rho) * (4. + 0.5 * (u**2-rho**2))
        a_2 /= (rho**2 * (4. + (u - rho)**2)**.5)
        a_3 = 2. * (u - rho)**2 * (1. + rho**2)
        a_3 /= (rho**2 * (u + rho) * (4. + (u - rho)**2)**.5)

        n = 4. * u * rho / (u + rho)**2
        k = 4. * n / (4. + (u - rho)**2)
        # We omit sqrt, because all python packages use k^2 convention.

        x_1 = ellipk(k)
        x_2 = ellipe(k)
        x_3 = ellip3(n, k)
        (x_1, x_2) = (x_2, x_1)  # WM94 under Eq. 9 are inconsitent with GR80.

        return (a_1*x_1 + a_2*x_2 + a_3*x_3) / np.pi

x = np.linspace(-3., 1., 100)
y = np.linspace(-3., 1., 100)

for x_ in x:
    for y_ in y:
        u = 10**x_
        rho = 10**y_
        out = _get_magnification_WM94(u, rho)
        print("{:e} {:e} {:e}".format(x_, y_, math.log10(out)))

