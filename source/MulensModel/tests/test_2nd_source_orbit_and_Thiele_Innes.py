#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing of 2nd source orbit calculation using Thiele Innes elements and transformation to
standard parametrization.
"""


import numpy as np
from numpy.testing import assert_almost_equal
from MulensModel.orbits.orbit import Orbit, OrbitEccentric
from PyAstronomy.pyasl import KeplerEllipse
from matplotlib import pyplot as plt
import MulensModel as mm


def on_line(A, B, C, eps=1e-3):
    x1, y1 = A
    x2, y2 = B
    x3, y3 = C
    pole = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
    if abs(pole) > eps:
        raise ValueError('periapsi error')
    return abs(pole) <= eps


def get_orbits(kwargs):
    """
    Return orbits in standard parametrization using both MulensModel and PyAstronomy
    """
    orbit_standard = OrbitEccentric(**kwargs)
    periapsis_epoch_standard = orbit_standard._periapsis_epoch
    orbit_PyAstronomy = KeplerEllipse(
        kwargs['semimajor_axis'],
        kwargs['period'],
        kwargs['eccentricity'],
        tau=periapsis_epoch_standard,
        w=kwargs['omega_periapsis'],
        i=kwargs['inclination'],
        Omega=kwargs['Omega_node'])
    return orbit_standard, orbit_PyAstronomy, periapsis_epoch_standard


def get_orbits_TI(kwargs):
    """
    Returns MulensModel orbit based on Thiele -Innes parametrization and trandsformation to standard parametrization
    """
    orbit_TI = Orbit(**kwargs)
    standard = orbit_TI.orbit_elements_dict_degrees()
    standard.pop('periapsis_epoch')
    return orbit_TI, standard


t_all_TI = 0.
t_all_Std = 0.

time_true = np.linspace(0, 3000, 1000)


# number of tests
n = 100


for i in range(n):
    # random generation of parameters
    t_0 = np.float64((np.random.rand() * (max(time_true)-min(time_true)))+min(time_true))
    parameters_TI = {
        't_0': t_0,
        'u_0': 0.0,
        't_E': 10000000000000000,
        'xi_period': 6000,
        'xi_A': np.float64(np.random.rand()*20)-10.,
        'xi_B': np.float64(np.random.rand()*20)-10.,
        'xi_F': np.float64(np.random.rand()*20)-10.,
        'xi_G': np.float64(np.random.rand()*20)-10.,
        'xi_argument_of_latitude_reference': np.float64(np.random.rand()*780)-360.,
        'xi_eccentricity': np.float64(np.random.rand()),
        'q_source': np.float64(np.random.rand())}

    zip_ = parameters_TI.items()
    orbit_parameters_TI = {key[3:]: value for (
        key, value) in zip_ if key[:3] == "xi_"}
    orbit_parameters_TI['epoch_reference'] = t_0

    orbit_TI, standard = get_orbits_TI(orbit_parameters_TI)

    # Primary source
    standard_1S = standard
    # Secondary source
    standard_2S = standard.copy()
    standard_2S['semimajor_axis'] /= parameters_TI['q_source']
    standard_2S['omega_periapsis'] += 180.
    standard_2S['argument_of_latitude_reference'] += 180.

    # Primary source
    orbit_standard_1, orbit_PyAstronomy_1, periapsis_epoch_standard_1 = get_orbits(
        standard_1S)
    # Primary orbit projection from MulensModel
    orbit_standard_1_xy = orbit_standard_1.get_reference_plane_position(
        time_true)
    # Primary orbit projection from PyAstronomy
    orbit_PyAstronomy_1_xy = orbit_PyAstronomy_1.xyzPos(time_true).T[:2]
    # Primary periapsis from MulensModel
    orbit_standard_1_xy_0 = orbit_standard_1.get_reference_plane_position(
        periapsis_epoch_standard_1)

    # Secondary source
    orbit_standard_2, orbit_PyAstronomy_2, periapsis_epoch_standard_2 = get_orbits(
        standard_2S)
    # Secondary orbit projection from MulensModel
    orbit_standard_2_xy = orbit_standard_2.get_reference_plane_position(
        time_true)
    # Secondary orbit projection from PyAstronomy
    orbit_PyAstronomy_2_xy = orbit_PyAstronomy_2.xyzPos(time_true).T[:2]
    # Secondary periapsis from MulensModel
    orbit_standard_2_xy_0 = orbit_standard_2.get_reference_plane_position(
        periapsis_epoch_standard_2)

    model_TI = mm.Model(parameters_TI)
    print('model:\n', parameters_TI)

    # shift from center of lens center of mass to center of mass of sources
    xi_delta = model_TI.parameters.source_1_parameters.xallarap_reference_position

    plt.scatter([0-xi_delta[0]], [0-xi_delta[1]], label='Sources CM TI',
                alpha=0.6, marker='.', s=25, c='red')

    # orbit from MulensModel based on standard parametrization
    plt.scatter(orbit_standard_1_xy[0]-xi_delta[0], orbit_standard_1_xy[1] -
                xi_delta[1], label='Std_MM_prime', alpha=0.6, marker='o', s=5, c='blue')
    plt.scatter(orbit_standard_2_xy[0]-xi_delta[0], orbit_standard_2_xy[1] -
                xi_delta[1], label='Std_MM_sec', alpha=0.6, marker='o', s=5, c='red')
    # periapsis from MulensModel based on standard parametrization
    plt.scatter(orbit_standard_1_xy_0[0]-xi_delta[0], orbit_standard_1_xy_0[1]-xi_delta[1],
                label='Std_MM_prime', alpha=0.6, marker='o', s=25, c='blue')
    plt.scatter(orbit_standard_2_xy_0[0]-xi_delta[0], orbit_standard_2_xy_0[1]-xi_delta[1],
                label='Std_MM_sec', alpha=0.6, marker='o', s=25, c='red')

    # orbit from PyAstronomy based on standard parametrization
    plt.scatter(orbit_PyAstronomy_1_xy[0]-xi_delta[0], orbit_PyAstronomy_1_xy[1] - xi_delta[1],
                label='PA_prime', alpha=0.6, marker="x", s=5)
    plt.scatter(orbit_PyAstronomy_2_xy[0]-xi_delta[0], orbit_PyAstronomy_2_xy[1] - xi_delta[1],
                label='PA_sec', alpha=0.6, marker="x", s=5)
    # orbits from  MulensModele based on Thiele-Innes parametrization
    plot_kwargs = {'t_start': min(time_true), 't_stop': max(time_true)}
    model_TI.plot_trajectory(label='both_sources TI', **plot_kwargs)
    orbit_TI_1_xy = orbit_TI.get_reference_plane_position(time_true)
    Periapsis_TI_1 = mm.Trajectory(
        times=[periapsis_epoch_standard_1], parameters=model_TI.parameters.source_1_parameters)
    Periapsis_TI_2 = mm.Trajectory(
        times=[periapsis_epoch_standard_2], parameters=model_TI.parameters.source_2_parameters)

    plt.scatter(Periapsis_TI_1.x, Periapsis_TI_1.y, label='TI_prime',
                alpha=0.6, marker='x', s=25, c='blue')
    plt.scatter(Periapsis_TI_2.x, Periapsis_TI_2.y, label='TI_sec',
                alpha=0.6, marker='x', s=25, c='red')
    plt.legend()
    plt.show()

    assert_almost_equal(orbit_standard_1_xy,
                        orbit_TI_1_xy, err_msg='Std_MM vs TI', decimal=2)
    assert_almost_equal(orbit_TI_1_xy, orbit_PyAstronomy_1_xy, err_msg='TI vs PA', decimal=2)
    assert_almost_equal(orbit_standard_1_xy, orbit_PyAstronomy_1_xy,
                        err_msg='RP1 vs PA1', decimal=2)
    assert_almost_equal(orbit_standard_2_xy, orbit_PyAstronomy_2_xy,
                        err_msg='RP2 vs PA2', decimal=2)
    print('CM, Periapsis_1 Periapsis_2 are aligned:', on_line(
        [Periapsis_TI_1.x, Periapsis_TI_1.y], [0-xi_delta[0], 0-xi_delta[1]], [Periapsis_TI_2.x, Periapsis_TI_2.y]))
