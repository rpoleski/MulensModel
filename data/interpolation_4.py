import math
import numpy as np
from math import sin, cos, sqrt
from scipy import integrate
from scipy.interpolate import interp1d, interp2d
from scipy.special import ellipk, ellipe
# These are complete elliptic integrals of the first and the second kind.
from sympy.functions.special.elliptic_integrals import elliptic_pi as ellip3

###################################################
# 1st and 2nd ellip integrals should be interpoalted between 1e-4 and 0.5
# ellip3(n, k) for n in (1e-3,0.4) and k as above 
###################################################

accuracy = 1.e-5
n_divide = 2 + 1
x_start = 1.e-4
x_stop = 0.5
y_start = 1.e-3
y_stop = 0.4
n_start = 10

def get_ellip(x, y):
    p = []
    for y_ in y:
        for x_ in x:
            index = (x_, y_)
            if index not in get_ellip.p:
                get_ellip.p[index] = ellip3(x_, y_)
            p.append(get_ellip.p[index])
    z = np.array(p).reshape(len(x), len(y))
    return z
get_ellip.p = dict()

x = np.logspace(np.log10(x_start), np.log10(x_stop), n_start)
y = np.logspace(np.log10(y_start), np.log10(y_stop), n_start)

iteration = 0 
add_x = [None]
add_y = [None]
while len(add_x) > 0 or len(add_y) > 0:
    iteration += 1
    add_x = []
    add_y = []
    p = get_ellip(x, y)

    interp_p = interp2d(x, y, p.T, kind='cubic')

    check_x = []
    for i in range(len(x)-1):
        check_x += np.logspace(np.log10(x[i]), np.log10(x[i+1]), n_divide)[1:-1].tolist()
    check_y = []
    for i in range(len(y)-1):
        check_y += np.logspace(np.log10(y[i]), np.log10(y[i+1]), n_divide)[1: -1].tolist()
    check_true_p = get_ellip(check_x, check_y).flatten()
    check_p = []
    check_x_ = []
    check_y_ = []
    for cx in check_x:
        for cy in check_y:
            check_p.append(interp_p(cx, cy)[0])
            check_x_.append(cx)
            check_y_.append(cy)
    check_p = np.array(check_p)
    relative_diff_p = np.abs(check_p - check_true_p) / check_true_p
    index = np.argsort(relative_diff_p)[-1]
    if relative_diff_p[index] < accuracy:
                continue
    add_x.append(check_x_[index])
    add_y.append(check_y_[index])
    new_x = np.sort(add_x + x.tolist())
    new_y = np.sort(add_y + y.tolist())
    print(iteration, len(new_x), len(new_y), np.max(relative_diff_p), add_x[0], add_y[0])
    x = new_x
    y = new_y

# TO DO:
# - check all log10 in interpolation_3.py
# - add log10 above
# - plot functions
