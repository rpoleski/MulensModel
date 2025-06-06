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


    parameters_extra = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt})

    return parameters_extra


def test_trajectory():
    """
    compares trajectory to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times(parameters)

    trajectory = mm.Trajectory(parameters=parameters, times=times)

    x_VBB = [-3.00059505, 1.59070972, 0.16689172, -1.25159717, -2.66308834]
    y_VBB = [-0.55177383, -0.16633580, -0.52931291, -0.99575576, -1.48861865]



    np.testing.assert_almost_equal(trajectory.x, x_VBB)
    np.testing.assert_almost_equal(trajectory.y, y_VBB)

