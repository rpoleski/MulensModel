"""
Calculates interpolation tables for elliptical integral of
the first and second kind.
"""
import math
import numpy as np
from math import sin, cos, sqrt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.special import ellipk, ellipe
# These are complete elliptic integrals of the first and the second kind.


accuracy = 1.e-6
n_divide = 10 + 1
x_start = 2.e-5
x_stop = 1. - 1.e-12
n_start = 10

# Settings end here.


def get_ellip(x):
    k = []
    e = []
    for x_ in x:
        if x_ not in get_ellip.k:
            get_ellip.k[x_] = ellipk(x_)
            get_ellip.e[x_] = ellipe(x_)
        k.append(get_ellip.k[x_])
        e.append(get_ellip.e[x_])
    return (np.array(k), np.array(e))
get_ellip.k = dict()
get_ellip.e = dict()

x = np.logspace(np.log10(x_start), np.log10(x_stop), n_start)

iteration = 0
add = [None]
while len(add) > 0:
    iteration += 1
    add = []
    (k, e) = get_ellip(x)
    interp_k = interp1d(np.log10(x), k, kind='cubic')
    interp_e = interp1d(np.log10(x), e, kind='cubic')
    for i in range(len(x)-1):
        x_1 = x[i]
        x_2 = x[i+1]
        check = np.logspace(np.log10(x_1), np.log10(x_2), n_divide)[1:-1]
        (check_true_k, check_true_e) = get_ellip(check)
        check_k = []
        check_e = []
        for c in check:
            check_k.append(interp_k(np.log10(c)))
            check_e.append(interp_e(np.log10(c)))
        check_k = np.array(check_k)
        check_e = np.array(check_e)
        relative_diff_k = np.abs(check_k - check_true_k) / check_true_k
        relative_diff_e = np.abs(check_e - check_true_e) / check_true_e
        relative_diff_max = np.amax(
            np.array([relative_diff_k, relative_diff_e]), axis=0)
        index = np.argsort(relative_diff_max)[-1]
        if relative_diff_max[index] < accuracy:
            continue
        add.append(check[index])
    new_x = np.sort(add + x.tolist())
    x = new_x

for (x_, k_, e_) in zip(x, k, e):
    print(x_, k_, e_)
