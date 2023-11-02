import unittest
import numpy as np

import MulensModel as mm


class TestModelGetMagnificationCurve(unittest.TestCase):
    """
    Test that Model can return magnification curves.
    """

    def setUp(self):
        self.parameters_source_1 = {'t_0': 5000., 'u_0': 0.005, 't_E': 25.}
        self.parameters_source_2 = {'t_0_2': 5100., 'u_0_2': 0.0003}
        self.times = np.arange(4900, 5200, 1.)
        self.mag_curve_source_1 = mm.MagnificationCurve(
            parameters=mm.ModelParameters(self.parameters_source_1),
            times=self.times)
        self.mag_curve_source_2 = mm.MagnificationCurve(
            parameters=mm.ModelParameters(self.parameters_source_2),
            times=self.times)

    def test_get_magnification_curve_1source(self):
        model = mm.Model(self.parameters_source_1)
        self.assertEqual(
            model.get_magnification_curve(self.times),
            self.mag_curve_source_1)

    def test_get_magnification_curve_2sources(self):
        parameters = self.parameters_source_2.copy()
        parameters['t_0_1'] = self.parameters_source_1['t_0']
        parameters['u_0_1'] = self.parameters_source_1['u_0']
        parameters['t_E'] = self.parameters_source_1['t_E']

        model = mm.Model(parameters)

        mag_curves_tuple_1 = model.get_magnification_curve(self.times)
        self.assertEqual(mag_curves_tuple_1[0], self.mag_curve_source_1)
        self.assertEqual(mag_curves_tuple_1[1], self.mag_curve_source_2)

        mag_curves_tuple_2 = model.get_magnification_curves(self.times)
        self.assertEqual(mag_curves_tuple_2[0], self.mag_curve_source_1)
        self.assertEqual(mag_curves_tuple_2[1], self.mag_curve_source_2)

        mag_curve_1 = model.get_magnification_curve_source_1(self.times)
        self.assertEqual(mag_curve_1, self.mag_curve_source_1)

        mag_curve_2 = model.get_magnification_curve_source_2(self.times)
        self.assertEqual(mag_curve_2, self.mag_curve_source_2)


class TestMagnificationCurve(unittest.TestCase):
    """
    Test that Magnification Curve stores more properties that can be accessed
    by the user.
    """

    def setUp(self):
        t_0 = 2456789.012345
        t_E = 23.4567
        u_0 = 1e-4
        rho = 1e-3
        self.gamma = 0.56789
        self.t_vec = np.array([-(rho ** 2 - u_0 ** 2) ** 0.5, 0.,
                          ((0.5 * rho) ** 2 - u_0 ** 2) ** 0.5])
        self.t_vec = self.t_vec * t_E + t_0

        self.params = mm.ModelParameters(
            {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})

        self.mag_curve = mm.MagnificationCurve(
            times=self.t_vec, parameters=self.params, gamma=self.gamma)
        methods = [t_0 - t_E, 'finite_source_uniform_Gould94', t_0 + t_E]
        self.mag_curve.set_magnification_methods(methods, 'point_source')

        self.pspl_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
        self.b0 = np.array([1.27323965, 0.19949906, 0.93421546])
        self.b1 = np.array([0.09489869, -0.03492121, -0.09655794])

    def test_A_pspl(self):
        self.assertEqual(
            self.mag_curve.A_pspl,
            self.pspl_model.get_magnification(self.t_vec))

    def test_B0B1(self):
        self.assertEqual(self.mag_curve.b0, self.b0)
        self.assertEqual(self.mag_curve.b1, self.b1)
        self.assertEqual(
            self.mag_curve.b0_gamma_b1, self.b0 - self.gamma * self.b1)

    def test_trajectory(self):
        traj = self.pspl_model.get_trajectory(self.t_vec)
        self.assertEqual(self.mag_curve.trajectory.x, traj.x)
        self.assertEqual(self.mag_curve.trajectory.y, traj.y)
