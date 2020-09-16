"""
Calculates interpolation tables for elliptical integral of the third kind.
"""
import math
import numpy as np
from math import sin, cos, sqrt
from scipy import integrate
from scipy.interpolate import interp1d, interp2d
from sympy.functions.special.elliptic_integrals import elliptic_pi as ellip3


accuracy = 1.e-4
n_divide = 3 + 1
x_start = 1.e-4
x_stop = 1.-1.e-6
y_start = 1.e-4
y_stop = 1.-1.e-6
n_start = 10
file_out_name = "interpolate_elliptic_integral_3.dat"

# Settings end here.


def get_ellip(x, y):
    p = []
    z = np.zeros((len(x), len(y)))
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

    interp_p = interp2d(x, y, p.T, kind='cubic')

    check_x = []
    for i in range(len(x)-1):
        check_ = np.logspace(np.log10(x[i]), np.log10(x[i+1]), n_divide)
        check_x += check_[1:-1].tolist()
    check_y = []
    for i in range(len(y)-1):
        check_ = np.logspace(np.log10(y[i]), np.log10(y[i+1]), n_divide)
        check_y += check_[1: -1].tolist()
    check_true_p = get_ellip(check_x, check_y)
    check_p = np.zeros((len(check_x), len(check_y)))
    for (ix, cx) in enumerate(check_x):
        for (iy, cy) in enumerate(check_y):
            if cy > cx:
                check_p[ix, iy] = 1.
                check_true_p[ix, iy] = 1.
            else:
                check_p[ix, iy] = interp_p(cx, cy)[0]
    relative_diff_p = np.abs(check_p - check_true_p) / check_true_p
    index = np.unravel_index(relative_diff_p.argmax(), relative_diff_p.shape)

    if np.max(relative_diff_p) < accuracy:
        continue
    add_x.append(check_x[index[0]])
    add_y.append(check_y[index[1]])
    new_x = np.sort(add_x + x.tolist())
    new_y = np.sort(add_y + y.tolist())
    x = new_x
    y = new_y

# Write to output file.
p = get_ellip(x, y)
with open(file_out_name, "w") as f_out:
    f_out.write(" ".join(["# X"] + [str(x_) for x_ in x] + ["\n"]))
    f_out.write(" ".join(["# Y"] + [str(y_) for y_ in y] + ["\n"]))
    for (i, x_) in enumerate(x):
        f_out.write(
            " ".join([str(p[i, j]) for (j, y_) in enumerate(y)] + ["\n"]))

# Read the output file and test it.
with open(file_out_name) as file_in:
    for line in file_in.readlines():
        if line[:3] == "# X":
            xx = np.array([float(t) for t in line.split()[2:]])
        if line[:3] == "# Y":
            yy = np.array([float(t) for t in line.split()[2:]])
pp = np.loadtxt(file_out_name)
print(np.all(x == xx))
print(np.all(y == yy))
print(np.all(p == pp))
