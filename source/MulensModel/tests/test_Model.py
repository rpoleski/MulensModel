import numpy as np
from numpy.testing import assert_almost_equal as almost
from math import isclose
import pytest
import unittest
import os.path

import MulensModel as mm

dir_2 = os.path.join(mm.DATA_PATH, 'unit_test_files')
dir_3 = os.path.join(mm.DATA_PATH, 'ephemeris_files')
SAMPLE_FILE_02_REF = os.path.join(dir_2, 'ob140939_OGLE_ref_v2.dat')  # HJD'
SAMPLE_FILE_03_EPH = os.path.join(dir_3, 'Spitzer_ephemeris_01.dat')  # UTC
SAMPLE_FILE_03_REF = os.path.join(dir_2, 'ob140939_Spitzer_ref_v2.dat')  # HJD'


def test_n_lenses():
    """check n_lenses property"""
    model_1 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.})
    model_2 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.,
                        "s": 1.1234, "q": 0.123, 'alpha': 192.34})
    model_3 = mm.Model({"t_0": 2456789., "u_0": 1., "t_E": 30.,
                        "s": 1.1234, "q": 0.123, 'alpha': 192.34,
                        'convergence_K': 0.04, 'shear_G': complex(0.1, -0.05)})
    assert model_1.n_lenses == 1
    assert model_2.n_lenses == 2
    assert model_3.n_lenses == 2


# Point Lens Tests
def test_model_PSPL_1():
    """tests basic evaluation of Paczynski model"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    times = np.array([t_0-2.5*t_E, t_0, t_0+t_E])
    model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    almost(model.get_magnification(times),
           np.array([1.028720763, 2.10290259, 1.26317278]),
           err_msg="PSPL model returns wrong values")


def test_model_init_1():
    """tests if basic parameters of Model.__init__() are properly passed"""
    t_0 = 5432.10987
    u_0 = 0.001
    t_E = 123.456
    rho = 0.0123
    my_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': rho})
    almost(my_model.parameters.t_0, t_0, err_msg='t_0 not set properly')
    almost(my_model.parameters.u_0, u_0, err_msg='u_0 not set properly')
    almost(my_model.parameters.t_E, t_E, err_msg='t_E not set properly')
    almost(my_model.parameters.rho, rho, err_msg='rho not set properly')


class TestModel(unittest.TestCase):
    def test_negative_t_E(self):
        with self.assertRaises(ValueError):
            mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': -100.})

    def test_finite_source_method_without_rho(self):
        """
        Test method that requires fintie source to run on model that
        has point source.
        """
        model = mm.Model({'t_0': 10, 'u_0': 0.2, 't_E': 50})
        with self.assertRaises(ValueError):
            model.set_magnification_methods(
                [9, "finite_source_uniform_Gould94", 11])

    def test_finite_source_method_without_rho_1(self):
        """
        2 source, only source 1 is finite, but we set finite method for both.
        """
        model = mm.Model({'t_0_1': 10, 'u_0_1': 0.2, 't_0_2': 30,
                          'u_0_2': 0.4, 't_E': 50, 'rho_1': 0.6})
        with self.assertRaises(ValueError):
            model.set_magnification_methods(
                [9, "finite_source_uniform_Gould94", 11])

    def test_finite_source_method_without_rho_2(self):
        """
        2 source, only source 1 is finite, but we set finite method for both.
        """
        model = mm.Model({'t_0_1': 10, 'u_0_1': 0.2, 't_0_2': 30,
                          'u_0_2': 0.4, 't_E': 50, 'rho_2': 0.6})
        with self.assertRaises(ValueError):
            model.set_magnification_methods(
                [29, "finite_source_uniform_Gould94", 31])


def test_model_methods():
    """
    Simplest model and setting point_source method.
    It obviously should work but it's worth to check it.
    """
    model = mm.Model({'t_0': 10, 'u_0': 0.2, 't_E': 50})
    model.set_magnification_methods([9, "point_source", 11])


class TestModelParallaxPIE(unittest.TestCase):
    def test_parallax_pi_E(self):
        """
        Make sure error is raised for API used in v1 and v2.
        """
        params = {'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'pi_E': (0.5, 0.6)}
        with self.assertRaises(KeyError):
            mm.ModelParameters(params)


def test_model_parallax_definition():
    """Update parameters in an existing model"""
    model_2 = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.,
                        'pi_E_N': 0.1, 'pi_E_E': 0.2})

    model_2.parameters.pi_E_N = 0.3
    model_2.parameters.pi_E_E = 0.4
    assert model_2.parameters.pi_E_N == 0.3
    assert model_2.parameters.pi_E_E == 0.4

    model_4 = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.,
                        'pi_E_N': 0.7, 'pi_E_E': 0.8})
    assert model_4.parameters.pi_E_N == 0.7
    assert model_4.parameters.pi_E_E == 0.8


def test_coords_transformation():
    """
    this was tested using http://ned.ipac.caltech.edu/forms/calculator.html
    """
    coords = "17:54:32.1 -30:12:34.0"
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.}, coords=coords)

    almost(model.coords.galactic_l.value, 359.90100049-360., decimal=4)
    almost(model.coords.galactic_b.value, -2.31694073, decimal=3)

    almost(model.coords.ecliptic_lon.value, 268.81102051, decimal=1)
    almost(model.coords.ecliptic_lat.value, -6.77579203, decimal=2)


def test_coords_input():
    """
    Test if wrong coords input raises an error.
    Only test that was missing to have 100% coverage in __init__.
    """
    ra, dec = "17:54:32.1", "-30:12:34.0"
    with pytest.raises(AttributeError):
        _ = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.}, ra=ra)
    with pytest.raises(AttributeError):
        _ = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100.}, dec=dec)


def test_init_parameters():
    """are parameters properly passed between Model and ModelParameters?"""
    t_0 = 6141.593
    u_0 = 0.5425
    t_E = 62.63
    params = mm.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    model = mm.Model(parameters=params)
    almost(model.parameters.t_0, t_0)
    almost(model.parameters.u_0, u_0)
    almost(model.parameters.t_E, t_E)


def test_limb_darkening():
    """check if limb_darkening coeffs are properly passed and converted"""
    gamma = 0.4555
    u = 3. * gamma / (2. + gamma)

    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'rho': 0.001})
    model.set_limb_coeff_gamma('I', gamma)

    almost(model.get_limb_coeff_gamma('I'), gamma)
    almost(model.get_limb_coeff_u('I'), u)


def test_limb_darkening_source():
    """check if limb_darkening is properly set using source=None"""
    gamma = 0.4555
    u = 3. * gamma / (2. + gamma)
    dict_1l2s = {'t_0_1': t_0, 'u_0_1': u_0, 't_0_2': t_0 + 50, 'u_0_2': u_0,
                 't_E': t_E, 'rho_1': 0.1, 'rho_2': 0.002}
    model = mm.Model(dict_1l2s)

    with pytest.raises(ValueError):
        model.get_limb_coeff_gamma('I', source=None)

    model.set_limb_coeff_gamma('I', gamma, source=None)
    assert model.get_limb_coeff_gamma('I', source=1) == gamma
    assert (model.get_limb_coeff_gamma('I') == [gamma]*model.n_sources)
    model.set_limb_coeff_u('I', u, source=1)
    assert model.get_limb_coeff_u('I', source=2) == u
    assert (model.get_limb_coeff_u('I') == [u]*model.n_sources)


def test_errors_in_get_magnification():
    """check if errors in get_magnification() are properly raised"""
    t_0 = 2450000.
    u_0 = 0.1
    dict_1l2s = {'t_0_1': t_0, 'u_0_1': u_0, 't_0_2': t_0 + 9, 'u_0_2': u_0, 't_E': 10., 'rho_1': 0.01, 'rho_2': 0.1}
    model = mm.Model(dict_1l2s)
    model.set_magnification_methods([2449999., "finite_source_uniform_Gould94", 2450001.], source=1)
    model.set_magnification_methods([2450008., "finite_source_uniform_Gould94", 2450010.], source=2)

    with pytest.raises(TypeError):
        model.get_magnification(2450000., source_flux_ratio=1)
    with pytest.raises(ValueError):
        model.get_magnification(2450000., separate=False)

    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'rho': 0.001})
    model.set_magnification_methods([2449999., "finite_source_uniform_Gould94", 2450001.])
    with pytest.raises(ValueError):
        model.get_magnification(2450000., source_flux_ratio=0.5)
    with pytest.raises(ValueError):
        model.get_magnification(2450000., separate=True)
    with pytest.raises(ValueError):
        model.get_magnification([np.nan])


def test_errors_in_limb_darkening():
    """check if limb_darkening errors are properly raised"""
    gamma = 0.4555
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'rho': 0.001})
    model.set_limb_coeff_gamma('I', gamma, source=None)

    with pytest.raises(ValueError):
        model.get_limb_coeff_gamma('I', source=3)
    with pytest.raises(ValueError):
        model.get_limb_coeff_u('I', source=3)

    times = np.arange(2449900., 2450101., 50)
    with pytest.raises(ValueError):
        model.get_magnification(times, bandpass='I', gamma=gamma)
    with pytest.raises(KeyError):
        model.get_magnification(times, bandpass='V')


def test_limb_darkening_with_uniform_methods_single_source():
    """
    Fail get_magnification() call with LD if no methods with LD are used,
    in case of a model with single source.
    """
    model = mm.Model({'t_0': 2450000., 'u_0': 0.1, 't_E': 100., 'rho': 0.001})
    model.set_limb_coeff_gamma('I', -1)
    model.default_magnification_method = 'finite_source_uniform_Gould94'
    times = np.arange(2449900., 2450101., 50)

    forbidden = {
        "finite_source_uniform_Gould94",
        "finite_source_uniform_Gould94_direct",
        "finite_source_uniform_WittMao94",
        "finite_source_uniform_Lee09",
    }
    for method in forbidden:
        model.set_magnification_methods([2449900., method, 2450100.])
        with pytest.raises(ValueError):
            model.get_magnification(times, bandpass='I')


def test_limb_darkening_with_uniform_methods_multi_source():
    """
    Fail get_magnification() with only uniform methods, for multiple source.
    LD and uniform method is set for one source only and the error persists.
    """
    t_0 = 2450000.
    u_0 = 0.1
    dict_1l2s = {'t_0_1': t_0, 'u_0_1': u_0, 't_0_2': t_0 + 50, 'u_0_2': u_0,
                 't_E': 100., 'rho_1': 0.001, 'rho_2': 0.1}
    model = mm.Model(dict_1l2s)
    model.set_limb_coeff_gamma('I', -1, source=1)
    model.default_magnification_method = 'finite_source_uniform_Gould94'

    times = np.arange(2449900., 2450101., 50)
    model.set_magnification_methods(
        [2449900., 'finite_source_uniform_WittMao94', 2450100.], source=1)
    with pytest.raises(ValueError):
        model.get_magnification(times, bandpass='I')

    model.set_limb_coeff_gamma('I', 2, source=2)
    model.set_magnification_methods(
        [2449950., 'finite_source_uniform_Lee09', 2450150.], source=2)
    with pytest.raises(ValueError):
        model.get_magnification(times, bandpass='I')


def test_different_limb_darkening():
    """testing effect of different limb darkening for different sources"""
    t_0 = 2450000.
    u_0 = 0.001
    t_E = 100.
    model_1 = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho': 0.1})
    model_2 = mm.Model({'t_0': t_0 + 50, 'u_0': u_0, 't_E': t_E, 'rho': 0.002})
    mag_method = [2449900., 'finite_source_LD_Yoo04', 2450100.]
    model_1.set_magnification_methods(mag_method)
    model_2.set_magnification_methods(mag_method)

    model_1.set_limb_coeff_gamma('I', -1)
    model_2.set_limb_coeff_gamma('I', 2)
    times = [2449990.0, 2449990.5, 2449991.0, 2450000.0, 2450010.0,
             2450048.5, 2450049.0, 2450049.9, 2450050.5, 2450051.0]
    mag_1 = model_1.get_magnification(times, bandpass='I')
    mag_2 = model_2.get_magnification(times, bandpass='I')

    dict_1l2s = {'t_0_1': t_0, 'u_0_1': u_0, 't_0_2': t_0 + 50, 'u_0_2': u_0,
                 't_E': t_E, 'rho_1': 0.1, 'rho_2': 0.002}
    model_1l2s = mm.Model(dict_1l2s)
    model_1l2s.set_limb_coeff_gamma('I', -1, source=1)
    model_1l2s.set_limb_coeff_gamma('I', 2, source=2)
    model_1l2s.set_magnification_methods(mag_method)
    mags = model_1l2s.get_magnification(times, bandpass='I', separate=True)

    assert (mag_1 == mags[0]).all()
    assert (mag_2 == mags[1]).all()


def test_set_magnification_methods_errors():
    """
    Testing if set_magnification_methods() works properly.
    This improved five lines of code coverage missing in model.py.
    """
    t_0 = 2450000.
    u_0 = 0.1
    dict_1l2s = {'t_0_1': t_0, 'u_0_1': u_0, 't_0_2': t_0 + 50, 'u_0_2': u_0,
                 't_E': 100., 'rho_1': 0.001, 'rho_2': 0.1}
    model = mm.Model(dict_1l2s)
    model_2 = mm.Model(dict_1l2s)
    methods = [2449900., 'point_source', 2450100.]

    with pytest.raises(TypeError):
        model.set_magnification_methods(2449900.)
    with pytest.raises(ValueError):
        model.set_magnification_methods(methods, source=0.1)
    with pytest.raises(ValueError):
        model.set_magnification_methods(methods, source=3)

    model.set_magnification_methods(methods, source=1)
    with pytest.raises(ValueError):
        model.set_magnification_methods(methods)
    model_2.set_magnification_methods(methods)
    with pytest.raises(ValueError):
        model_2.set_magnification_methods(methods, source=1)


def test_t_E():
    """make sure t_E can be accessed properly"""
    t_0 = 2460000.
    u_0 = 0.1
    t_E = 12.3456
    t_star = 0.01234
    params = dict(t_0=t_0, u_0=u_0, t_E=t_E)
    model_1 = mm.Model(params)
    params['t_star'] = t_star
    model_2 = mm.Model(params)

    almost(model_1.parameters.t_E, t_E)
    almost(model_2.parameters.t_E, t_E)


# Binary Lens tests
# Binary lens parameters:
alpha = 49.58
s = 1.3500
q = 0.00578
# Other parameters
t_E = 62.63
rho = 0.01

# Adjust t_0, u_0 from primary to CM
shift_x = -s * q / (1. + q)
delta_t_0 = -t_E * shift_x * np.cos(alpha * np.pi / 180. + np.pi)
delta_u_0 = -shift_x * np.sin(alpha * np.pi / 180. + np.pi)
t_0 = 2456141.593 + delta_t_0
u_0 = 0.5425 + delta_u_0


def test_BLPS_01():
    """simple binary lens with point source"""
    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q})

    model = mm.Model(parameters=params)
    t = np.array([2456112.5])
    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    magnification = model.get_magnification(data.time[0])
    almost(magnification[0], 4.691830781584699)
# This value comes from early version of this code.
# almost(m, 4.710563917)
# This value comes from Andy's getbinp().


class TestBLPS_ShearActive(unittest.TestCase):
    """
    simple binary lens with point source and external convergence and shear
    """
    def setUp(self):
        self.params = mm.ModelParameters({
            't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
            'q': q, 'convergence_K': 0.08, 'shear_G': complex(0.1, -0.1)})

        self.model = mm.Model(parameters=self.params)
        self.t = np.array([2456112.5])
        self.data = mm.MulensData(data_list=[self.t, self.t*0.+16.,
                                             self.t*0.+0.01])

    def test_vbbl_fail(self):
        self.model.default_magnification_method = 'point_source'
        magnification = self.model.get_magnification(self.data.time[0])
        assert not isclose(magnification[0], 4.691830781584699, abs_tol=1e-2)

    def test_wm95_fail(self):
        self.model.default_magnification_method = 'point_source_WM95'
        magnification = self.model.get_magnification(self.data.time[0])
        assert not isclose(magnification[0], 4.691830781584699, abs_tol=1e-2)


class TestBLPS_Shear(unittest.TestCase):
    """
    simple binary lens with point source and external convergence and shear
    """

    def setUp(self):
        self.params = mm.ModelParameters({
            't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
            'q': q, 'convergence_K': 0.0, 'shear_G': complex(0, 0)})

        self.model = mm.Model(parameters=self.params)
        self.t = np.array([2456112.5])
        self.data = mm.MulensData(data_list=[self.t, self.t*0.+16.,
                                             self.t*0.+0.01])

    def test_vbbl(self):
        self.model.default_magnification_method = 'point_source'
        magnification = self.model.get_magnification(self.data.time[0])
        almost(magnification[0], 4.691830781584699)

    def test_wm95(self):
        self.model.default_magnification_method = 'point_source_WM95'
        magnification = self.model.get_magnification(self.data.time[0])
        almost(magnification[0], 4.691830781584699)


def test_BLPS_02():
    """
    simple binary lens with extended source and different methods to
    evaluate magnification
    """

    params = mm.ModelParameters({
        't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
        'q': q, 'rho': rho})
    model = mm.Model(parameters=params)

    t = np.array([6112.5, 6113., 6114., 6115., 6116., 6117., 6118., 6119])
    t += 2450000.
    methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole', 2456116.5,
               'VBBL', 2456117.5]
    model.set_magnification_methods(methods)
    assert model.default_magnification_method == 'point_source'
    assert model.methods == methods

    data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])
    result = model.get_magnification(data.time)

    expected = np.array([4.69183078, 2.87659723, 1.83733975, 1.63865704,
                         1.61038135, 1.63603122, 1.69045492, 1.77012807])
    almost(result, expected, decimal=4)

    # Possibly, this test should be re-created in test_FitData.py
    # Below we test passing the limb coeff to VBBL function.
    # data.bandpass = 'I'
    model.set_limb_coeff_u('I', 10.)
    # This is an absurd value but I needed something quick.
    result = model.get_magnification(
        data.time, gamma=model.get_limb_coeff_gamma('I'))
    almost(result[5], 1.6366862, decimal=3)
    result_2 = model.get_magnification(data.time, bandpass='I')
    almost(result, result_2)


class TestBLPS02AC(unittest.TestCase):
    """
    simple binary lens with extended source and different methods to evaluate
    magnification - version with adaptivecontouring
    """

    def setUp(self):
        t = np.array([6112.5, 6113., 6114., 6115., 6116., 6117., 6118., 6119])
        t += 2450000.
        self.data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])

        ac_name = 'Adaptive_Contouring'
        methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                   2456116.5,
                   ac_name, 2456117.5]
        accuracy_1 = {'accuracy': 0.04}
        accuracy_2 = {'accuracy': 0.01, 'ld_accuracy': 0.00001}

        params = mm.ModelParameters({
            't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
            'q': q, 'rho': rho})
        self.model_ac_1 = mm.Model(parameters=params)
        self.model_ac_1.set_magnification_methods(methods)
        self.model_ac_1.set_magnification_methods_parameters(
            {ac_name: accuracy_1})

        self.model_ac_2 = mm.Model(parameters=params)
        # data.bandpass = 'I'
        self.model_ac_2.set_limb_coeff_u('I', 10.)
        # This is an absurd value but I needed something quick.
        self.model_ac_2.set_magnification_methods(methods)
        self.model_ac_2.set_magnification_methods_parameters(
            {ac_name: accuracy_2})

    def test_methods(self):
        def test_model_methods(model):
            assert model.default_magnification_method == 'point_source'
            methods_compare = [
                2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                2456116.5, 'Adaptive_Contouring', 2456117.5]
            assert model.methods == methods_compare

        test_model_methods(self.model_ac_1)
        test_model_methods(self.model_ac_2)

    def test_methods_parameters_1(self):
        # test get_magnification_methods_parameters()
        assert (self.model_ac_1.get_magnification_methods_parameters(
            'Adaptive_Contouring') ==
                {'adaptive_contouring': {'accuracy': 0.04}})

    def test_mag_calculation_1(self):
        """Test calculation of magnification"""
        result = self.model_ac_1.get_magnification(self.data.time)
        expected = np.array([4.69183078, 2.87659723, 1.83733975, 1.63865704,
                             1.61038135, 1.63603122, 1.69045492, 1.77012807])
        almost(result, expected, decimal=3)

    def test_methods_parameters_2(self):
        """
        test get_magnification_methods_parameters()
        and methods_parameters()
        """
        AC = 'Adaptive_Contouring'
        dict_1 = self.model_ac_2.get_magnification_methods_parameters(AC)
        reference = {AC.lower(): {'accuracy': 0.01, 'ld_accuracy': 0.00001}}
        assert dict_1 == reference

    def test_mag_calculation_2(self):
        """Test calculation of magnification with limb darkening"""
        result = self.model_ac_2.get_magnification(
            self.data.time, gamma=self.model_ac_2.get_limb_coeff_gamma('I'))
        almost(result[5], 1.6366862, decimal=3)

    def test_wrong_parameter_name(self):
        """
        Test if KeyError is raised when wrong parameter name is used.
        """
        ac_name = 'Adaptive_Contouring'
        wrong_par = {ac_name: {'acuracy': 0.04}}
        with self.assertRaises(KeyError):
            self.model_ac_2.set_magnification_methods_parameters(wrong_par)
        wrong_par_2 = {ac_name: {'accuracy': 0.01, 'ld_acuracy': 1e-5}}
        with self.assertRaises(KeyError):
            self.model_ac_2.set_magnification_methods_parameters(wrong_par_2)


class TestMethodsParameters(unittest.TestCase):
    """
    make sure additional parameters are properly passed to very inner functions
    """
    def setUp(self):
        t = np.array([2456117.])
        self.data = mm.MulensData(data_list=[t, t*0.+16., t*0.+0.01])

        self.params = mm.ModelParameters({
            't_0': t_0, 'u_0': u_0, 't_E': t_E, 'alpha': alpha, 's': s,
            'q': q, 'rho': rho})
        methods = [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                   2456116.5, 'VBBL', 2456117.5]
        self.model_1 = mm.Model(parameters=self.params)
        self.model_1.set_magnification_methods(methods)

        vbbl_options_2 = {'accuracy': 0.1}
        methods_parameters_2 = {'VBBL': vbbl_options_2}
        self.model_2 = mm.Model(parameters=self.params)
        self.model_2.set_magnification_methods(methods)
        self.model_2.set_magnification_methods_parameters(methods_parameters_2)

        vbbl_options_3 = {'accuracy': 1.e-5}
        methods_parameters_3 = {'VBBL': vbbl_options_3}
        self.model_3 = mm.Model(parameters=self.params)
        self.model_3.set_magnification_methods(methods)
        self.model_3.set_magnification_methods_parameters(methods_parameters_3)

    def test_methods(self):
        def test_model_methods(model):
            assert (model.default_magnification_method == 'point_source')
            assert (model.methods ==
                    [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                     2456116.5,
                     'VBBL', 2456117.5])

        test_model_methods(self.model_1)
        test_model_methods(self.model_2)
        test_model_methods(self.model_3)

    def test_get_magnification_methods(self):
        assert (self.model_1.get_magnification_methods() ==
                [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                 2456116.5,
                 'VBBL', 2456117.5])
        assert (self.model_1.get_magnification_methods(source=1) ==
                [2456113.5, 'Quadrupole', 2456114.5, 'Hexadecapole',
                 2456116.5,
                 'VBBL', 2456117.5])
        with self.assertRaises(IndexError):
            self.model_1.get_magnification_methods(source=2)

    def test_mag_calculations(self):
        result_1 = self.model_1.get_magnification(self.data.time)
        result_2 = self.model_2.get_magnification(self.data.time)
        result_3 = self.model_3.get_magnification(self.data.time)

        assert result_1[0] != result_2[0]
        assert result_1[0] != result_3[0]
        assert result_2[0] != result_3[0]

    def test_get_magnification_methods_parameters(self):
        with self.assertRaises(KeyError):
            self.model_1.get_magnification_methods_parameters('vbbl')

        assert (self.model_2.get_magnification_methods_parameters(
            'vbbl') == {'vbbl': {'accuracy': 0.1}})
        assert (self.model_3.get_magnification_methods_parameters(
            'vbbl') == {'vbbl': {'accuracy': 1.e-5}})

    def test_default_magnification_methods(self):
        """
        Test if methods are properly changed, no deprecated method.
        """
        model = mm.Model(self.params)
        assert model.default_magnification_method == 'point_source'

        model.default_magnification_method = 'point_source_point_lens'
        assert model.default_magnification_method == 'point_source_point_lens'

        model.default_magnification_method = 'VBBL'
        assert model.default_magnification_method == 'VBBL'

    def test_wrong_parameter_name(self):
        """
        Test if KeyError is raised when wrong parameter name is used.
        """
        wrong_par = {'VBBL': {'acuracy': 1.e-5}}
        with self.assertRaises(KeyError):
            self.model_3.set_magnification_methods_parameters(wrong_par)


def test_caustic_for_orbital_motion():
    """
    check if caustics calculated for different epochs in orbital motion model
    are as expected
    """
    q = 0.1
    s = 1.3
    params = {'t_0': 100., 'u_0': 0.1, 't_E': 10., 'q': q,
              's': s, 'ds_dt': 0.5, 'alpha': 0., 'dalpha_dt': 0.}
    model = mm.Model(parameters=params)

    model.update_caustics()
    almost(model.caustics.get_caustics(),
           mm.CausticsBinary(q=q, s=s).get_caustics())

    model.update_caustics(100.+365.25/2)
    almost(model.caustics.get_caustics(),
           mm.CausticsBinary(q=q, s=1.55).get_caustics())


def test_update_single_lens_with_shear_caustic():
    """
    make sure that updating single lens caustic works ok
    """
    convergence_K = 0.1
    shear_G = complex(-0.1, -0.2)

    model = mm.Model(mm.ModelParameters({
        't_0': 0., 'u_0': 1., 't_E': 2., 'alpha': 183.,
        'convergence_K': 0., 'shear_G': complex(0, 0)}))
    model.parameters.convergence_K = convergence_K
    model.parameters.shear_G = shear_G
    model.update_caustics()
    assert model.caustics.convergence_K == convergence_K
    assert model.caustics.shear_G == shear_G


def test_magnifications_for_orbital_motion():
    """
    make sure that orbital motion parameters are properly passed to
    magnification methods calculations
    """
    dict_static = {'t_0': 100., 'u_0': 0.1, 't_E': 100., 'q': 0.99,
                   's': 1.1, 'alpha': 190.}
    dict_motion = dict_static.copy()
    dict_motion.update({'ds_dt': -2, 'dalpha_dt': -300.})
    static = mm.Model(dict_static)
    motion = mm.Model(dict_motion)

    t_1 = 100.
    almost(static.get_magnification(t_1), motion.get_magnification(t_1))

    t_2 = 130.
    static.parameters.s = 0.93572895
    static.parameters.alpha = 165.359342916
    almost(static.get_magnification(t_2), motion.get_magnification(t_2))


def test_model_binary_and_finite_sources():
    """
    test if model magnification calculation for binary source works with
    finite sources (both rho and t_star given)
    """
    model = mm.Model({
        't_0_1': 5000., 'u_0_1': 0.005, 'rho_1': 0.001,
        't_0_2': 5100., 'u_0_2': 0.0003, 't_star_2': 0.03, 't_E': 25.})
    model_1 = mm.Model(model.parameters.source_1_parameters)
    model_2 = mm.Model(model.parameters.source_2_parameters)

    (t1, t2) = (4999.95, 5000.05)
    (t3, t4) = (5099.95, 5100.05)
    finite = 'finite_source_uniform_Gould94'
    methods = [t1, finite, t2, 'point_source', t3, finite, t4]
    model.set_magnification_methods(methods)
    model_1.set_magnification_methods(methods)
    model_2.set_magnification_methods(methods)

    def test_model_methods(test_model):
        assert (test_model.default_magnification_method == 'point_source')
        assert (test_model.methods ==
                [t1, finite, t2, 'point_source', t3, finite, t4])

    test_model_methods(model)
    test_model_methods(model_1)
    test_model_methods(model_2)

    (f_s_1, f_s_2) = (100., 300.)
    time = np.linspace(4900., 5200., 4200)
    mag_1 = model_1.get_magnification(time)
    mag_2 = model_2.get_magnification(time)

    # test: model.set_source_flux_ratio(f_s_2/f_s_1)
    fitted = model.get_magnification(time, source_flux_ratio=f_s_2 / f_s_1)
    expected = (mag_1 * f_s_1 + mag_2 * f_s_2) / (f_s_1 + f_s_2)
    almost(fitted, expected)

    # test separate=True option:
    (mag_1_, mag_2_) = model.get_magnification(time, separate=True)
    almost(mag_1, mag_1_)
    almost(mag_2, mag_2_)


def test_binary_source_and_fluxes_for_bands():
    """
    Test if setting different flux ratios for different bands in binary
    source models works properly. The argument flux_ratio_constraint
    is set as string.
    """
    model = mm.Model({'t_0_1': 5000., 'u_0_1': 0.05,
                      't_0_2': 5100., 'u_0_2': 0.003, 't_E': 30.})

    times_I = np.linspace(4900., 5200., 3000)
    times_V = np.linspace(4800., 5300., 250)
    (f_s_1_I, f_s_2_I) = (10., 20.)
    (f_s_1_V, f_s_2_V) = (15., 5.)
    q_f_I = f_s_2_I / f_s_1_I
    q_f_V = f_s_2_V / f_s_1_V
    (mag_1_I, mag_2_I) = model.get_magnification(times_I, separate=True)
    (mag_1_V, mag_2_V) = model.get_magnification(times_V, separate=True)
    effective_mag_I = (mag_1_I + mag_2_I * q_f_I) / (1. + q_f_I)
    effective_mag_V = (mag_1_V + mag_2_V * q_f_V) / (1. + q_f_V)
    # flux_I = mag_1_I * f_s_1_I + mag_2_I * f_s_2_I + f_b_I
    # flux_V = mag_1_V * f_s_1_V + mag_2_V * f_s_2_V + f_b_V

    # model.set_source_flux_ratio_for_band('I', q_f_I)
    # model.set_source_flux_ratio_for_band('V', q_f_V)

    # Test Model.get_magnification()
    result_I = model.get_magnification(times_I, source_flux_ratio=q_f_I)
    result_V = model.get_magnification(times_V, source_flux_ratio=q_f_V)
    almost(result_I, effective_mag_I)
    almost(result_V, effective_mag_V)


class TestSeparateMethodForEachSource(unittest.TestCase):
    """
    Checks if setting separate magnification method for each source in
    binary source models works properly.
    """
    def setUp(self):
        parameters = {'t_0_1': 5000., 'u_0_1': 0.01, 'rho_1': 0.005,
                      't_0_2': 5100., 'u_0_2': 0.001, 'rho_2': 0.005,
                      't_E': 1000.}
        self.model = mm.Model(parameters)

    def test_1(self):
        self.model.set_magnification_methods(
            [4999., 'finite_source_uniform_Gould94', 5101.], source=1)
        # In order not to get "no FS method" warning:
        self.model.set_magnification_methods(
            [0., 'finite_source_uniform_Gould94', 1.], source=2)
        out = self.model.get_magnification(5000., separate=True)
        almost([out[0][0], out[1][0]], [103.46704167, 10.03696291])
        assert self.model.default_magnification_method == 'point_source'
        assert (self.model.methods[1] ==
                [4999., 'finite_source_uniform_Gould94', 5101.])
        assert (self.model.methods[2] ==
                [0., 'finite_source_uniform_Gould94', 1.])
        assert (self.model.get_magnification_methods(source=1) ==
                [4999., 'finite_source_uniform_Gould94', 5101.])
        assert (self.model.get_magnification_methods(source=2) ==
                [0., 'finite_source_uniform_Gould94', 1.])
        assert (self.model.get_magnification_methods() ==
                {1: [4999., 'finite_source_uniform_Gould94', 5101.],
                 2: [0., 'finite_source_uniform_Gould94', 1.]})

    def test_2(self):
        self.model.set_magnification_methods(
            [4999., 'finite_source_uniform_Gould94', 5001.], source=1)
        self.model.set_magnification_methods(
            [5099., 'finite_source_uniform_Gould94', 5101.], source=2)
        out = self.model.get_magnification(5100., separate=True)
        almost([out[0][0], out[1][0]], [9.98801936, 395.96963727])
        assert self.model.default_magnification_method == 'point_source'
        assert (self.model.methods[1] ==
                [4999., 'finite_source_uniform_Gould94', 5001.])
        assert (self.model.methods[2] ==
                [5099., 'finite_source_uniform_Gould94', 5101.])


def test_get_lc():
    """
    Test if Model.get_lc() works properly; we test on binary source model
    without finite source effect.
    """
    model = mm.Model({'t_0_1': 5000., 'u_0_1': 1.,
                      't_0_2': 5100., 'u_0_2': 0.1,
                      't_E': 100.})
    out = model.get_lc(5050., source_flux=[1., 2.], blend_flux=3.)
    almost(out, 19.668370500043526)


def test_is_finite_source():
    model_fs = mm.Model({'t_0': 10, 'u_0': 1, 't_E': 3, 'rho': 0.001})
    model_ps = mm.Model({'t_0': 10, 'u_0': 1, 't_E': 3})

    assert model_fs.parameters.is_finite_source()
    assert not model_ps.parameters.is_finite_source()


def test_repr():
    """Test if printing is Model instance is OK."""
    parameters = {'t_0': 2454656.4, 'u_0': 0.003,
                  't_E': 11.1, 't_star': 0.055}
    begin = ("    t_0 (HJD)       u_0    t_E (d)    t_star (d) \n"
             "2454656.40000  0.003000    11.1000      0.055000 \n")
    end = "default magnification method: point_source"
    model = mm.Model(parameters)
    assert str(model) == begin + end

    coords = "17:54:32.10 -30:12:34.99"
    model = mm.Model(parameters, coords=coords)
    expected = "{:}coords: {:}\n{:}".format(begin, coords, end)
    assert str(model) == expected

    model = mm.Model(parameters)
    methods = [2454656.3, 'finite_source_uniform_Gould94', 2454656.5]
    model.set_magnification_methods(methods)
    expected = "{:}{:}\nother magnification methods: {:}".format(
        begin, end, methods)
    assert str(model) == expected

    model = mm.Model(parameters)
    model.set_limb_coeff_gamma("I", 0.5)
    expected = begin + end + "\nlimb-darkening coeffs (gamma): [{'I': 0.5}]"
    assert str(model) == expected


def prepare_xallarap_test(
        xi_Omega_node=0., xi_argument_of_latitude_reference=0.,
        t_0_xi=None, xi_eccentricity=None, xi_omega_periapsis=None):
    """
    prepare data for unit tests of xallarap models
    """
    t_0 = 2459876.0
    d_time = 4.
    u_0 = 0.01357
    xi_a = 0.0123456789
    tau = 0.03838383

    common = {'t_0': t_0, 't_E': d_time / tau, 'u_0': u_0}
    xallarap = {
        'xi_period': 2.*d_time, 'xi_semimajor_axis': xi_a,
        'xi_Omega_node': xi_Omega_node, 'xi_inclination': 0.,
        'xi_argument_of_latitude_reference': xi_argument_of_latitude_reference
        }
    if t_0_xi is not None:
        if t_0_xi == "t_0":
            t_0_xi = t_0
        elif t_0_xi == "t_0+d_time":
            t_0_xi = t_0 + d_time
        xallarap['t_0_xi'] = t_0_xi
    if xi_eccentricity is not None:
        xallarap['xi_eccentricity'] = xi_eccentricity
    if xi_omega_periapsis is not None:
        xallarap['xi_omega_periapsis'] = xi_omega_periapsis

    model_1 = mm.Model(common)
    model_2 = mm.Model({**common, **xallarap})
    return (model_1, model_2, t_0, d_time, tau, u_0, xi_a)


def test_xallarap_at_t_0_circular():
    """
    Make sure that xallarap circular and non-xallarap 1L1S models
    produce the same magnifications at t_0.
    """
    (model_1, model_2, t_0) = prepare_xallarap_test()[:2+1]

    assert model_1.get_magnification(t_0) == model_2.get_magnification(t_0)


def test_xallarap_at_t_0_eccentric():
    """
    Make sure that xallarap eccentric and non-xallarap 1L1S models
    produce the same magnifications at t_0.
    """
    kwargs = {'xi_eccentricity': 0.5, 'xi_omega_periapsis': 12.3456}
    (model_1, model_2, t_0) = prepare_xallarap_test(**kwargs)[:2+1]

    assert model_1.get_magnification(t_0) == model_2.get_magnification(t_0)


def test_xallarap_at_t_0_plus_half_of_period_1():
    """
    Xallarap - circular orbit, half period after t_0, and Omega+nu_0 = 0.
    Expected u is from pen and pencil calculations.
    """
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test()[1:]

    u2 = u_0**2 + (tau - 2. * xi_a)**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_2():
    """
    Xallarap - circular orbit, half period after t_0, and Omega+nu_0 = 90.
    Expected u is from pen and pencil calculations.
    """
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(90., 0.)[1:]

    u2 = (u_0 - 2. * xi_a)**2 + tau**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_3():
    """
    Xallarap - circular orbit, half period after t_0, and Omega+nu_0 = 180.
    Expected u is from pen and pencil calculations.
    """
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(90., 90.)[1:]

    u2 = u_0**2 + (tau + 2. * xi_a)**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_4():
    """
    Xallarap - circular orbit, half period after t_0, and Omega+nu_0 = 180.
    The t_0_xi is provided as a parameter and = t_0.
    Expected u is from pen and pencil calculations.
    """
    args = [90., 90., "t_0"]
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(*args)[1:]

    u2 = u_0**2 + (tau + 2. * xi_a)**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_5():
    """
    Xallarap - circular orbit, half period after t_0, and Omega+nu_0 = 180.
    The t_0_xi is provided as a parameter and = t_0.
    Expected u is from pen and pencil calculations.
    """
    args = [90., 90., "t_0+d_time"]
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(*args)[1:]

    u2 = u_0**2 + tau**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_6_eccentric():
    """
    Extremely eccentric xallarap orbit checked half period later
    """
    kwargs = {'xi_eccentricity': 0.999, 'xi_omega_periapsis': 0.}
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(**kwargs)[1:]

    u2 = u_0**2 + (tau - 2. * xi_a)**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    assert expected == model.get_magnification(t_0+d_time)


def test_xallarap_at_t_0_plus_half_of_period_7_eccentric():
    """
    Extremely eccentric xallarap orbit with Omega=270 checked half period later
    """
    kwargs = {'xi_eccentricity': 0.999, 'xi_omega_periapsis': 0.,
              'xi_Omega_node': 270.}
    (model, t_0, d_time, tau, u_0, xi_a) = prepare_xallarap_test(**kwargs)[1:]

    u2 = (u_0 + 2. * xi_a)**2 + tau**2
    expected = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
    almost(expected, model.get_magnification(t_0+d_time))


class TestGetTrajectory(unittest.TestCase):

    def setUp(self):
        # Parallax model parameters
        self.model_parameters_par = {
            't_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
            'pi_E_N': -0.248, 'pi_E_E': 0.234, 't_0_par': 2456836.2}
        self.coords = "17:47:12.25 -21:22:58.2"

        self.model_with_par = mm.Model(
            self.model_parameters_par, coords=self.coords)
        self.model_with_par.parallax(satellite=True, earth_orbital=True,
                                     topocentric=False)

        self.ref_OGLE = np.loadtxt(SAMPLE_FILE_02_REF, unpack=True)
        self.times_OGLE = self.ref_OGLE[0] + 2450000.
        self.ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True)
        self.ephemerides_file = SAMPLE_FILE_03_EPH
        self.times_spz = self.ref_Spitzer[0] + 2450000.

    def test_1L1S(self):
        # straight-up trajectory for static point-lens model
        t_0 = self.model_parameters_par['t_0']
        u_0 = self.model_parameters_par['u_0']
        t_E = self.model_parameters_par['t_E']
        times = np.arange(t_0 - t_E, t_0 + t_E, 1.)

        model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
        trajectory = model.get_trajectory(times)
        y = np.ones(len(times)) * u_0
        x = (times - t_0) / t_E
        ratio_x = trajectory.x / x
        ratio_y = trajectory.y / y
        np.testing.assert_almost_equal(ratio_x, 1.)
        np.testing.assert_almost_equal(ratio_y, 1.)

    def test_1L1S_annual_parallax(self):
        """case with annual parallax"""
        trajectory = self.model_with_par.get_trajectory(
            self.times_OGLE)

        ratio_x = trajectory.x / self.ref_OGLE[6]
        ratio_y = trajectory.y / self.ref_OGLE[7]
        np.testing.assert_almost_equal(ratio_x, [1.] * len(ratio_x), decimal=3)
        np.testing.assert_almost_equal(ratio_y, [1.] * len(ratio_y), decimal=3)

    def test_1L1S_satellite_parallax_1(self):
        """Case with satellite parallax (check test_Model_Parallax.py)"""
        satellite_skycoord_obj = mm.SatelliteSkyCoord(
            ephemerides_file=self.ephemerides_file)
        satellite_skycoord = satellite_skycoord_obj.get_satellite_coords(
            self.times_spz)
        trajectory = self.model_with_par.get_trajectory(
            self.times_spz, satellite_skycoord=satellite_skycoord)

        ratio_x = trajectory.x / self.ref_Spitzer[6]
        ratio_y = trajectory.y / self.ref_Spitzer[7]
        np.testing.assert_almost_equal(ratio_x, [1.] * len(ratio_x), decimal=2)
        np.testing.assert_almost_equal(ratio_y, [1.] * len(ratio_y), decimal=3)

    def test_1L1S_satellite_parallax_2(self):
        """Case with satellite parallax (check test_Model_Parallax.py)"""
        model_with_sat_par = mm.Model(
            self.model_parameters_par, ephemerides_file=self.ephemerides_file,
            coords=self.coords)
        trajectory = model_with_sat_par.get_trajectory(self.times_spz)

        ratio_x = trajectory.x / self.ref_Spitzer[6]
        ratio_y = trajectory.y / self.ref_Spitzer[7]
        np.testing.assert_almost_equal(ratio_x, [1.] * len(ratio_x), decimal=2)
        np.testing.assert_almost_equal(ratio_y, [1.] * len(ratio_y), decimal=3)

    def test_1L2S(self):
        """Binary source trajectories"""
        model = mm.Model({
            't_0_1': 5000., 'u_0_1': 0.005, 'rho_1': 0.001,
            't_0_2': 5100., 'u_0_2': 0.0003, 't_star_2': 0.03, 't_E': 25.})
        model_1 = mm.Model(model.parameters.source_1_parameters)
        model_2 = mm.Model(model.parameters.source_2_parameters)

        time = np.linspace(4900., 5200., 4200)

        (traj_1_1L2S, traj_2_1L2S) = model.get_trajectory(time)
        traj_1_1L1S = model_1.get_trajectory(time)
        np.testing.assert_equal(traj_1_1L2S.x, traj_1_1L1S.x)
        np.testing.assert_equal(traj_1_1L2S.y, traj_1_1L1S.y)

        traj_2_1L1S = model_2.get_trajectory(time)
        np.testing.assert_equal(traj_2_1L2S.x, traj_2_1L1S.x)
        np.testing.assert_equal(traj_2_1L2S.y, traj_2_1L1S.y)


class TestNSources(unittest.TestCase):

    def setUp(self):
        """
        Parameters correspond to OB151459: Hwang et al. 2018, AJ, 155, 259
        Data are simulated without noise.
        Fluxes in paper were in the 18th magnitude system. MM works with
        MAG_ZEROPOINT=22, so these magnitudes have an offset relative to the
        original event.
        Table 3: Static 1L3S Model
        """
        self.t_E = 4.921

        self.source_1_params = {'t_E': self.t_E, 't_0': 7199.946, 'u_0': 0.065}
        self.source_2_params = {'t_E': self.t_E, 't_0': 7200.193, 'u_0': 2.638e-3, 'rho': 4.503e-3}
        self.source_3_params = {'t_E': self.t_E, 't_0': 7200.202, 'u_0': 0.281e-3, 'rho': 0.631e-3}

        self.flux_ratios = {'q_2': 0.014, 'q_3': 0.006}
        self.source_flux_1 = 0.056
        self.source_flux_2 = self.source_flux_1 * self.flux_ratios['q_2']
        self.source_flux_3 = self.source_flux_1 * self.flux_ratios['q_3']
        self.source_fluxes = [self.source_flux_1, self.source_flux_2, self.source_flux_3]
        self.blend_flux = 0.100

        self.mag_methods = [7200.1, 'finite_source_uniform_Gould94', 7200.3]

        self.model_1 = mm.Model(self.source_1_params)
        self.model_2 = mm.Model(self.source_2_params)
        self.model_2.set_magnification_methods(self.mag_methods)
        self.model_3 = mm.Model(self.source_3_params)
        self.model_3.set_magnification_methods(self.mag_methods)
        self.models = [self.model_1, self.model_2, self.model_3]

        self.times = [7199.946, 7200., 7200.193, 7200.202]

        self.magnifications = np.vstack((
            self.model_1.get_magnification(self.times),
            self.model_2.get_magnification(self.times),
            self.model_3.get_magnification(self.times)))

        self.flux = self.source_flux_1 * self.model_1.get_magnification(self.times)
        self.flux += self.source_flux_2 * self.model_2.get_magnification(self.times)
        self.flux += self.source_flux_3 * self.model_3.get_magnification(self.times)
        self.flux += self.blend_flux

        self.model_params = {'t_E': self.t_E}
        for i, source_params in enumerate(
                [self.source_1_params, self.source_2_params, self.source_3_params]):
            for key, value in source_params.items():
                if key != 't_E':
                    self.model_params['{0}_{1}'.format(key, i+1)] = value

        self.model = mm.Model(self.model_params)
        self.model.set_magnification_methods(self.mag_methods, 2)
        self.model.set_magnification_methods(self.mag_methods, 3)

    def test_get_magnification(self):
        np.testing.assert_almost_equal(self.model.get_magnification(self.times), self.magnifications, decimal=4)

    def test_get_lc_1(self):
        """
        Returns :
            magnification: *numpy.ndarray*
                Magnification values for each epoch.
        """
        model_mags = self.model.get_lc(self.times, source_flux=self.source_fluxes, blend_flux=self.blend_flux)
        np.testing.assert_almost_equal(model_mags, mm.Utils.get_mag_from_flux(self.flux), decimal=4)

    def test_get_lc_2(self):
        """
        Returns :
            magnification: *numpy.ndarray*
                Magnification values for each epoch.
        """
        model_mags = self.model.get_lc(
            self.times, source_flux=self.source_flux_1,
            source_flux_ratio=[self.flux_ratios['q_2'], self.flux_ratios['q_3']], blend_flux=self.blend_flux)
        np.testing.assert_almost_equal(
            model_mags, mm.Utils.get_mag_from_flux(self.flux), decimal=4
        )

    def test_get_magnification_curves(self):
        mag_curves = self.model.get_magnification_curves(self.times, None, None)
        for mag_curve, model in zip(mag_curves, self.models):
            np.testing.assert_almost_equal(
                mag_curve.get_magnification(), model.get_magnification(self.times), decimal=6
            )

    def test_set_times(self):
        times = self.model.set_times()
        np.testing.assert_almost_equal(
            self.source_1_params['t_0'] - 1.5 * self.t_E, times[0], decimal=4
        )
        np.testing.assert_almost_equal(
            self.source_3_params['t_0'] + 1.5 * self.t_E, times[-1], decimal=4
        )

# def test_N_sources_gamma():
#    """
#    Test a model with gammas for different sources.
#    """
#    raise NotImplementedError()


# Tests to Add:
#
# test set_times:
#   keywords to test:
#     t_range=None, t_start=None, t_stop=None, dt=None, n_epochs=1000
#
# test set_default_magnification_method:
#   change from default value
#
# test get_satellite_coords: (check test_Model_Parallax.py)
#   returns None if no ephemerides file set
#   other condidtions probably covered by other unit tests
#
# test _magnification_2_sources: check instances of q_flux being specified vs.
# separate. Specifically worried that calls to magnification from other parts
# of the code work as expected.
#
# properties: parallax, caustics, parameters, n_lenses, n_source, is_static,
# coords, bandpasses,
