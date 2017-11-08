import unittest
import numpy as np

from astropy import units as u

from MulensModel.modelparameters import ModelParameters


class TestModelParameters(unittest.TestCase):
    def test_too_many_parameters_for_init(self):
        with self.assertRaises(ValueError):
            mp = ModelParameters(
                {'pi_E':(1., 1.), 'pi_E_N':1.})
        with self.assertRaises(ValueError):
            mp = ModelParameters(
                {'pi_E':(1., 1.), 'pi_E_E':1.})

def test_init_parameters():
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_E':t_E})

    np.testing.assert_almost_equal(params.t_0, t_0)
    np.testing.assert_almost_equal(params.u_0, u_0)
    np.testing.assert_almost_equal(params.t_E, t_E.value)

def test_repr_parameters():
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63*u.day
    params = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_E':t_E})
    
    out_1 = "  t_0 (HJD)       u_0    t_E (d) \n"
    out_2 = " 6141.59300  0.542500    62.6300 \n"
    
    assert (out_1 + out_2) == str(params)

def test_rho_t_e_t_star():
    """check if conversions between rho, t_E, and t_star work ok"""
    t_0 = 2450000.
    u_0 = 0.1
    t_E = 20. * u.day
    rho = 0.001
    t_star = t_E * rho

    params_1 = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_E':t_E, 'rho':rho})
    np.testing.assert_almost_equal(params_1.t_star.value, t_star.value)

    params_2 = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_star':t_star, 'rho':rho})
    np.testing.assert_almost_equal(params_2.t_E.value, t_E.value)

    params_3 = ModelParameters({'t_0':t_0, 'u_0':u_0, 't_star':t_star, 't_E':t_E})
    np.testing.assert_almost_equal(params_3.rho.value, rho)

class test(unittest.TestCase):
    def test_too_much_rho_t_e_t_star(self):
        with self.assertRaises(ValueError):
            t_0 = 2450000.
            u_0 = 0.1
            t_E = 20. * u.day
            rho = 0.001
            t_star = t_E * rho
            ModelParameters({'t_0':t_0, 'u_0':u_0, 't_E':t_E, 'rho':rho, 't_star':t_star})

