import numpy as np
import os
import unittest

import MulensModel as mm


SAMPLE_FILE = os.path.join(mm.DATA_PATH, 'unit_test_files', 'FSPL_test_1.dat')

class PointLensTestCase(unittest.TestCase):

    def setUp(self):
        with open(SAMPLE_FILE) as data_file:
            lines = data_file.readlines()
            ulens_params = lines[2].split()

        self.parameters = mm.ModelParameters(
                {'t_0': float(ulens_params[1]), 'u_0': float(ulens_params[2]),
                 't_E': float(ulens_params[3]), 'rho': float(ulens_params[4])})
        self.gamma = float(ulens_params[5])

        names = ['Time', 'b_0', 'b_1', 'Mag_FS', 'Mag_LD', 'Mag']
        self.data = np.genfromtxt(SAMPLE_FILE, names=names)

        tau = (self.data['Time'] - self.parameters.t_0) / self.parameters.t_E
        self.u = np.sqrt(self.parameters.u_0**2 + tau**2)
        self.z = self.u / self.parameters.rho
        self.pspl_magnification = (self.u**2 + 2.) / (
                self.u * np.sqrt(self.u**2 + 4.))

        self.point_lens_obj = mm.PointLens(parameters=self.parameters)

    def test_B_0_function(self):
        """test private _B_0_function"""
        test_b_0 = self.point_lens_obj._B_0_function(self.z)
        np.testing.assert_almost_equal(test_b_0, self.data['b_0'], decimal=5)

    def test_B_1_function(self):
        """test private _B_1_function"""
        test_b_1 = self.point_lens_obj._B_1_function(self.z)
        np.testing.assert_almost_equal(test_b_1, self.data['b_1'], decimal=4)

    def test_get_point_lens_finite_source_magnification(self):
        """test PLFS"""
        test_FSPL = self.point_lens_obj.get_point_lens_finite_source_magnification(
            self.u, self.pspl_magnification)
        np.testing.assert_almost_equal(test_FSPL, self.data['Mag_FS'], decimal=5)

    def test_get_point_lens_limb_darkening_magnification(self):
        """test PLFS+LD"""
        test_FSPL_LD = self.point_lens_obj.get_point_lens_limb_darkening_magnification(
            self.u, self.pspl_magnification, self.gamma)
        np.testing.assert_almost_equal(test_FSPL_LD / self.data['Mag_LD'], 1.,
                                       decimal=4)
