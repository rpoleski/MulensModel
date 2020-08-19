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
n_divide = 5 + 1
x_start = 1.e-4
x_stop = 0.5
y_start = 1.e-3
y_stop = 0.4
n_start = 10

def get_ellip(x, y):
    p = []
    z = np.zeros( (len(x), len(y)) )
    for (i, x_) in enumerate(x):
        for (j, y_) in enumerate(y):
            index = (x_, y_)
            if index not in get_ellip.p:
                get_ellip.p[index] = ellip3(x_, y_)
            z[i, j] = get_ellip.p[index]
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

    #interp_p = interp2d(x, y, p.T, kind='cubic')
    interp_p = interp2d(np.log10(x), np.log10(y), p.T, kind='cubic')

    check_x = []
    for i in range(len(x)-1):
        check_x += np.logspace(np.log10(x[i]), np.log10(x[i+1]), n_divide)[1:-1].tolist()
    check_y = []
    for i in range(len(y)-1):
        check_y += np.logspace(np.log10(y[i]), np.log10(y[i+1]), n_divide)[1: -1].tolist()
    check_true_p = get_ellip(check_x, check_y)
    check_p = np.zeros( (len(check_x), len(check_y)) )
    for (ix, cx) in enumerate(check_x):
        for (iy, cy) in enumerate(check_y):
            #check_p[ix, iy] = interp_p(cx, cy)[0]
            check_p[ix, iy] = interp_p(np.log10(cx), np.log10(cy))[0]
    relative_diff_p = np.abs(check_p - check_true_p) / check_true_p
    index = np.unravel_index(relative_diff_p.argmax(), relative_diff_p.shape)

    if np.max(relative_diff_p) < accuracy:
                continue
    add_x.append(check_x[index[0]])
    add_y.append(check_y[index[1]])
    new_x = np.sort(add_x + x.tolist())
    new_y = np.sort(add_y + y.tolist())
    print(iteration, len(new_x), len(new_y), np.max(relative_diff_p), add_x[0], add_y[0])
    x = new_x
    y = new_y

p = get_ellip(x, y)
for (i, x_) in enumerate(x):
    for (j, y_) in enumerate(y):
        #print(x_, y_, p[i, j])
        pass

