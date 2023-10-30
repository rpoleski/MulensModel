import numpy as np
import os
import unittest

import MulensModel as mm


SAMPLE_FILE = os.path.join(mm.DATA_PATH, 'unit_test_files', 'FSPL_test_1.dat')


def calc_u_and_z(parameters, time):
    tau = (time - parameters.t_0) / parameters.t_E
    u_ = np.sqrt(parameters.u_0 ** 2 + tau ** 2)
    z_ = u_ / parameters.rho
    return (u_, z_)


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

        (self.u, self.z) = calc_u_and_z(self.parameters, self.data['Time'])
        self.pspl_magnification = (self.u**2 + 2.) / (
                self.u * np.sqrt(self.u**2 + 4.))

        self.point_lens_obj = mm.PointLens(parameters=self.parameters)
        self.FS_uniform_obj = mm.PointLensFiniteSourceUniformGould94(
            parameters=self.parameters)
        self.FS_LD_obj = mm.PointLensFiniteSourceLDYoo04(
            parameters=self.parameters, gamma=self.gamma)

    def test_B_0_function(self):
        """test private _B_0_function"""
        for obj in [self.point_lens_obj, self.FS_uniform_obj, self. FS_LD_obj]:
            test_b_0 = obj._B_0_function(self.z)
            np.testing.assert_almost_equal(test_b_0, self.data['b_0'], decimal=5)

    def test_B_1_function(self):
        """test private _B_1_function"""
        for obj in [self.point_lens_obj, self.FS_LD_obj]:
            test_b_1 = obj._B_1_function(self.z)
            np.testing.assert_almost_equal(test_b_1, self.data['b_1'], decimal=5)

        with self.assertRaises(NameError):
            test_b_1 = self.FS_uniform_obj._B_1_function(self.z)

    def test_get_point_lens_point_source_magnification(self):
        for obj in [self.point_lens_obj, self.FS_uniform_obj, self. FS_LD_obj]:
            test_pspl_mag = obj.get_point_lens_point_source_magnification()
            np.testing.assert_almost_equal(
                test_pspl_mag, self.pspl_magnification, decimal=5)
            test_pspl_mag_1 = obj.get_pspl_magnification()
            np.testing.assert_almost_equal(
                test_pspl_mag_1, self.pspl_magnification, decimal=5)

    def test_get_point_lens_finite_source_magnification(self):
        """test PLFS"""
        # rename:
        # point_lens_finite_source
        # -->
        # finite_source_uniform_Gould94
        test_FSPL_0 = self.point_lens_obj.get_finite_source_uniform_Gould94_magnification(
            self.u)
        np.testing.assert_almost_equal(
            test_FSPL_0, self.data['Mag_FS'], decimal=5)

        # Which one(s)?
        test_FSPL = self.FS_uniform_obj.get_fspl_magnification(
            self.u)
        # test_FSPL = obj.get_finite_source_magnification(
        #     self.u)
        # test_FSPL = obj.get_finite_source_point_lens_magnification(
        #     self.u)
        np.testing.assert_almost_equal(
            test_FSPL, self.data['Mag_FS'], decimal=5)

    def test_get_point_lens_limb_darkening_magnification(self):
        """test PLFS+LD"""
        # rename:
        # point_lens_limb_darkening
        # -->
        # finite_source_LD_Yoo94
        test_FSPL_0 = self.point_lens_obj.get_finite_source_LD_Yoo94_magnification(
            self.u, self.gamma)
        np.testing.assert_almost_equal(
            test_FSPL_0 / self.data['Mag_LD'], 1., decimal=4)

        test_FSPL_LD = self.FS_LD_obj.get_fspl_magnification(self.u)
        np.testing.assert_almost_equal(
            test_FSPL_LD / self.data['Mag_LD'], 1., decimal=4)


# Adapted from test_MagnificationCurve
class Test_Lee09_and_WittMao94(unittest.TestCase):
    """
    test Lee et al. 2009 and Witt & Mao 1994 finite source calculation
    """

    def setUp(self):
        self.t_vec = np.array([3.5, 2., 1., 0.5, 0.])

        self.params_0 = mm.ModelParameters(
            {'t_0': 0., 'u_0': 0.5, 't_E': 1., 'rho': 1.})
        self.params_1 = mm.ModelParameters(
            {'t_0': 0., 'u_0': 0.1, 't_E': 1., 'rho': 1.})

        (self.u_0, _) = calc_u_and_z(self.params_0, self.t_vec)
        (self.u_1, _) = calc_u_and_z(self.params_1, self.t_vec)

        self.point_lens_0 = mm.PointLens(self.params_0)
        self.point_lens_1 = mm.PointLens(self.params_1)

        # The values below were calculated using code developed by P. Mroz.
        self.expected_0 = np.array(
            [1.01084060513, 1.06962639343, 1.42451408166, 2.02334097551,
             2.13919086656])
        self.expected_1 = np.array(
            [1.01110609638, 1.07461016241, 1.57232954942, 2.21990790526,
             2.39458814753])
        #    expected_2 = np.array([1.0110829794, 1.07404148634, 1.55620547462,
        #                           2.24809136704, 2.44503143812])
        # The last values are for 2-parameter LD with same settings and lambda=0.3.
        # Correction is:
        #  -lambda*(1-1.25*sqrt(costh))
        # and for 1-parameter LD we used:
        #  1-gamma*(1-1.5*costh)

    def test_uniform_source(self):
        # Test uniform source first.
        test_0 = mm.PointLensFiniteSourceUniformLee09(
            parameters=self.params_0)
        results_0 = test_0.get_fspl_magnification(self.u_0)
        np.testing.assert_almost_equal(self.expected_0, results_0, decimal=4)

        # rename:
        # point_lens_uniform_integrated
        # -->
        # finite_source_uniform_Lee09
        results_2 = self.point_lens_0.get_finite_source_uniform_Lee09_magnification(
            self.u_0)
        np.testing.assert_almost_equal(self.expected_0, results_2, decimal=4)

    def test_1parm_LD(self):
        # Then test 1-parameter limb-darkening.
        test_1 = mm.PointLensFiniteSourceLDLee09(
            parameters=self.params_1)
        results_1 = test_1.get_fspl_magnification(self.u_1)
        np.testing.assert_almost_equal(self.expected_1, results_1, decimal=3)

        # rename:
        # point_lens_LD_integrated
        # -->
        # finite_source_LD_Lee09
        results_2 = self.point_lens_1.get_finite_source_LD_Lee09_magnification(
            self.u_1, self.gamma)
        np.testing.assert_almost_equal(self.expected_1, results_2, decimal=4)

    def test_WittMao94_0(self):
        # Tests for Witt & Mao 1994 start here
        test_0 = mm.PointLensFiniteSourceUniformWittMao94(
            parameters=self.params_0)
        results_0 = test_0.get_fspl_magnification(self.u_0)
        np.testing.assert_almost_equal(self.expected_0, results_0, decimal=4)

        # rename:
        # point_lens_large_finite_source
        # -->
        # finite_source_uniform_WittMao94
        results_2 = self.point_lens_0.get_finite_source_uniform_WittMao94_magnification(
            self.u_0)
        np.testing.assert_almost_equal(self.expected_0, results_2, decimal=4)

    def test_WittMao94_1(self):
        test_1 = mm.PointLensFiniteSourceLDWittMao94(
            parameters=self.params_1)
        results_3 = test_1.get_fspl_magnification(self.u_1)
        np.testing.assert_almost_equal(self.expected_1, results_3, decimal=3)

        # rename:
        # point_lens_large_LD_integrated
        # -->
        # finite_source_LD_WittMao94
        results_2 = self.point_lens_1.get_finite_source_LD_WittMao94_magnification(
            self.u_1, self.gamma)
        np.testing.assert_almost_equal(self.expected_1, results_2, decimal=4)
