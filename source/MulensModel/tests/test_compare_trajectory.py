import MulensModel as mm
import numpy as np
import unittest
import pytest

def get_times(parameters):
    N = 5
    times = []
    t_0 = parameters.t_0
    t_E = parameters.t_E

    for i in range(0,N,1):
        times.append(t_0 - 3 * t_E + i * (6 * t_E / (N - 1)))
    return times

def get_parameters():

    s = 1.2
    q = 0.123
    u_0 = 0.555
    alpha = 17.5
    t_0 = 2455500
    t_E = 100
    rho = 0.01

    ds_dt = 20.2
    dalpha_dt = -30.3
    ds_z_dt = 20
    a = 2


    parameters_extra = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt})

    return parameters_extra


def test_separation_for_circular_orbit():
    """
    compares separation to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times(parameters)

    model = mm.Model(parameters=parameters)

    separation = model.parameters.get_s(times)

    separation_VBB_circular = [0.57896196, 0.36839104, 1.20000000, 1.66168431, 1.61032808]

    np.testing.assert_almost_equal(separation, separation_VBB_circular)

def test_trajectory_for_circular_orbit():
    """
    compares trajectory to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times(parameters)

    trajectory = mm.Trajectory(parameters=parameters, times=times)

    x_VBB_circular = [-3.00057308, 1.59071791, 0.16689172, -1.25159957, -2.66309601]
    y_VBB_circular = [-0.55189331, -0.16625747, -0.52931291, -0.99575273, -1.48860493]

    np.testing.assert_almost_equal(trajectory.x, x_VBB_circular)
    np.testing.assert_almost_equal(trajectory.y, y_VBB_circular)

def test_separation_for_elliptical_orbit():
    """
    compares separation to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times(parameters)

    model = mm.Model(parameters=parameters)

    separation = model.parameters.get_s(times)

    separation_VBB_elliptical = [0.64550215, 1.42405587, 1.20000000, 1.34447172, 2.14780495]

    np.testing.assert_almost_equal(separation, separation_VBB_elliptical)

def test_trajectory_for_elliptical_orbit():
    """
    compares trajectory to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times(parameters)

    trajectory = mm.Trajectory(parameters=parameters, times=times)

    x_VBB_elliptical = [-3.03964927, 1.59911138, 0.16689172, 1.23132793, 2.67836603]
    y_VBB_elliptical = [-0.26183446, -0.02945825, -0.52931291, 1.02071373, 1.46095188]


    np.testing.assert_almost_equal(trajectory.x, x_VBB_elliptical)
    np.testing.assert_almost_equal(trajectory.y, y_VBB_elliptical)

