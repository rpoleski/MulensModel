import unittest
import numpy as np

from astropy import units as u

from MulensModel.modelparameters import ModelParameters


class TestModelParameters(unittest.TestCase):
    def test_too_many_parameters_for_init(self):
        with self.assertRaises(ValueError):
            mp = ModelParameters(pi_E=(1., 1.), pi_E_N=1.)
        with self.assertRaises(ValueError):
            mp = ModelParameters(pi_E=(1., 1.), pi_E_E=1.)

def test_init_parameters():
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = ModelParameters(t_0=t_0, u_0=u_0, t_E=t_E)

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E.value)

