import unittest

import MulensModel as mm


def test_single_lens_with_mass_sheet():
    """
    Test if Chang-Refsdal microlensing parameters are properly defined.
    """
    basic = {'t_0': 1000., 'u_0': 0.1, 't_E': 20.}
    G = complex(-0.1, -0.2)
    K = -0.1

    parameters = mm.ModelParameters({**basic})
    parameters = mm.ModelParameters({**basic, 'shear_G': G, 'alpha': 123.})
    parameters = mm.ModelParameters({**basic, 'convergence_K': K})
    parameters = mm.ModelParameters(
        {**basic, 'shear_G': G, 'convergence_K': K, 'alpha': 123.})


class TestParameters(unittest.TestCase):
    def test_failing_single_lens_with_mass_sheet(self):
        """
        Test if Chang-Refsdal microlensing fails when it's expected to fail.
        """
        basic = {'t_0': 1000., 'u_0': 0.1, 't_E': 20.}
        G = complex(-0.1, -0.2)
        K = -0.1
        alpha = 123.

        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters(
                {'t_0': 1000., 'u_0': 0.1, 'shear_G': G})
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters(
                {'t_E': 20., 'u_0': 0.1, 'shear_G': G})
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters(
                {'t_0': 1000., 't_E': 20., 'shear_G': G})

        # Cases below have too many parameters:
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'alpha': alpha})
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters(
                {**basic, 'convergence_K': K, 'alpha': alpha})
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({
                **basic, 'convergence_K': K, 'alpha': alpha,
                'dalpha_dt': -0.3})
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({
                **basic, 'shear_G': G, 'convergence_K': K, 'alpha': alpha,
                'dalpha_dt': -0.3})

        # The case below is missing alpha:
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'shear_G': G})

