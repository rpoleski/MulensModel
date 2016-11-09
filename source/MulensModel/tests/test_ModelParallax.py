import unittest

from MulensModel.modelparameters import ModelParameters


class TestModelParallax(unittest.TestCase):
    def test_too_many_parameters_for_init(self):
        with self.assertRaises(ValueError):
            mp = ModelParameters(pi_E=(1., 1.), pi_E_N=1.)
        with self.assertRaises(ValueError):
            mp = ModelParameters(pi_E=(1., 1.), pi_E_E=1.)

