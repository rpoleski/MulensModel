import MulensModel as mm
import numpy as np
import unittest
import pytest


def test_separation():
    """
    compares trajectory to values form VBBinaryLensing
    """


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

    N = 5
    times = []

    for i in range(0,N,1):
        times.append(t_0 - 3 * t_E + i * (6 * t_E / (N - 1)))

    parameters_extra = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha,'rho': rho, 'ds_dt': ds_dt, 'dalpha_dt': dalpha_dt, 'ds_z_dt': ds_z_dt})

    trajectory = mm.Trajectory(parameters=parameters_extra, times=times)

    model = mm.Model(parameters=parameters_extra)
  
    separation = model.parameters.get_s(times)
    
    separation_VBB = [0.918651, 1.332512, 1.200000, 1.079128, 1.427853]
    np.testing.assert_almost_equal(separation, separation_VBB)

 
