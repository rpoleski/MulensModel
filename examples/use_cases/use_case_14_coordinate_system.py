"""
Use Case 14 - jcy

MulensModel will be written assuming the coordinate system
(t0,u0,alpha) are defined relative to the center of mass. This is not
always the most efficient choice for fitting. This use case covers the
conversion from center of magnification coordinates to center of mass
coordinates. Essentially, it is the responsibility of the user to
convert from their MCMC coordinate system to the center of mass
coordinate system needed for the magnification calculation.

"""
import numpy as np

import MulensModel as mm


raise NotImplementedError('frame_origin not implemented for Model')


def convert_cof_mag2mass(t0, te, u0, alpha, s, q):
    """
    function to convert from center of magnification to center of mass
    coordinates. Note that this function is for illustration only. It has
    not been tested and may have sign errors.
    """
    if s <= 1.0:
        return t0, u0
    else:
        delta = q / (1. + q) / s
        delta_u0 = delta * np.sin(alpha * np.pi / 180.)
        delta_tau = delta * np.cos(alpha * np.pi / 180.)
        t0_prime = t0 + delta_tau * te
        u0_prime = u0 + delta_u0
        return t0_prime, u0_prime


# Define model parameters in CoMAGN system
t0_center_of_mag = 7000.
u0_center_of_mag = 0.1
alpha_center_of_mag = 30.
te = 30.

print('Center of magnification: {0}, {1}'.format(
        t0_center_of_mag, u0_center_of_mag))

s = 1.1
q = 0.001

# Get parameters in CoMASS system
(t0_center_of_mass, u0_center_of_mass) = convert_cof_mag2mass(
        t0_center_of_mag, te, u0_center_of_mag, alpha_center_of_mag, s, q)

print('Center of mass: {0}, {1}'.format(t0_center_of_mass, u0_center_of_mass))

# How does this get passed to a minimizer?

# Alternatively,
model = mm.Model(
            {'t_0': 2457000., 'u_0': 0.1, 't_E': 30., 'rho': 0.001,
             'alpha': 30, 's': 1.1, 'q': 0.001},
            frame_origin='magnification')

print(model.parameters.t_0, model.parameters.u_0)
