import unittest

import MulensModel as mm


t0 = 1000
u0 = 0.1
tE = 20
K = -0.1
G = complex(-0.1, -0.2)
alpha = 123.
da_dt = -0.3

basic = {'t_0': t0, 'u_0': u0, 't_E': tE}

def test_1():
    parameters = mm.ModelParameters({**basic})

def test_2():
    parameters = mm.ModelParameters({**basic, 'shear_G': G})

def test_3():
    parameters = mm.ModelParameters({**basic, 'convergence_K': K, 'alpha': alpha})

def test_4():
    parameters = mm.ModelParameters({**basic, 'shear_G': G, 'convergence_K': K, 'alpha': alpha})

# In each case below, calling mm.ModelParameters() should raise an error.
class TestParameters(unittest.TestCase):
    def test_5(self):
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({'t_0': t0, 'u_0': u0, 'shear_G': G})

    def test_6(self):
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'alpha': alpha})

    def test_7(self):
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'shear_G': G, 'alpha': alpha})

    def test_8(self):
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'convergence_K': K, 'alpha': alpha, 'dalpha_dt': da_dt})

    def test_9(self):
        with self.assertRaises(KeyError):
            parameters = mm.ModelParameters({**basic, 'shear_G': G, 'convergence_K': K, 'alpha': alpha, 'dalpha_dt': da_dt})

