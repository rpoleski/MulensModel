import numpy as np
import unittest

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


def configure():
    """
    Prepare PointLensWithShear and Trajectory instances for tests
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

    return (point_lens, trajectory)


def test_get_ps_with_shear_magnification_2():
    """
    Test magnification calculation for point lens with convergence.
    """
    (point_lens, trajectory) = configure()
    test_pspl_shear = point_lens.get_point_source_magnification(trajectory)
    np.testing.assert_almost_equal(test_pspl_shear[0], 11.7608836, decimal=5)


class test_errors(unittest.TestCase):
    def test_function_input(self):
        """
        Make sure input of the main function is Trajectory type
        """
        (lens, trajectory) = configure()
        with self.assertRaises(TypeError):
            lens.get_point_source_magnification(123.)

    def test_not_implemented_functions(self):
        """
        Make sure that finite source functions inherited from PointLens
        are not implemented.
        """
        (lens, trajectory) = configure()

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_finite_source_magnification(trajectory)

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_limb_darkening_magnification(trajectory)

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_uniform_integrated_magnification(trajectory)

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_LD_integrated_magnification(trajectory)

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_large_finite_source_magnification(trajectory)

        with self.assertRaises(NotImplementedError):
            lens.get_point_lens_large_LD_integrated_magnification(trajectory)
