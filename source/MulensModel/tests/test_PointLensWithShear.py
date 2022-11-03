import numpy as np
import os

import MulensModel as mm


def test_get_ps_with_shear_magnification_1():
    """test PLPS+KG"""
    t_0 = 1000.
    t_E = 20.
    u_0 = 0.1
    t_vec = np.array([0., 100.]) * t_E + t_0
    convergence_K = 0.1
    shear_G = complex(-0.1, 0.2)
    parameters = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E,
        'convergence_K': convergence_K, 'shear_G': shear_G,
        'alpha': 0.})
    point_lens = mm.PointLensWithShear(parameters=parameters)
    # Set trajectory to be a single point
    trajectory = mm.Trajectory(t_vec, parameters)
    test_pspl_shear = point_lens.get_point_source_magnification(trajectory)
    np.testing.assert_almost_equal(test_pspl_shear[0], 5.556327, decimal=5)


def test_get_ps_with_shear_magnification_2():
    """
    Test magnification calculation for point lens with convergence.
    """
    t_0 = 1000.
    t_E = 20.
    u_0 = 0.1
    t_vec = np.array([0., 100.]) * t_E + t_0
    convergence_K = 0.1

    parameters = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'convergence_K': convergence_K})
    point_lens = mm.PointLensWithShear(parameters=parameters)
    trajectory = mm.Trajectory(t_vec, parameters)
    test_pspl_shear = point_lens.get_point_source_magnification(trajectory)
    np.testing.assert_almost_equal(test_pspl_shear[0], 11.7608836, decimal=5)
