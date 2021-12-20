import os
import numpy as np
import unittest
from astropy import units as u
from astropy.coordinates import SkyCoord

import MulensModel as mm


dir_1 = os.path.join(mm.DATA_PATH, 'photometry_files', 'OB140939')
dir_2 = os.path.join(mm.DATA_PATH, 'unit_test_files')
dir_3 = os.path.join(mm.DATA_PATH, 'ephemeris_files')

SAMPLE_FILE_02 = os.path.join(dir_1, 'ob140939_OGLE.dat')  # HJD'
SAMPLE_FILE_02_REF = os.path.join(dir_2, 'ob140939_OGLE_ref_v1.dat')  # HJD'
SAMPLE_FILE_03 = os.path.join(dir_1, 'ob140939_Spitzer.dat')  # HJD'
SAMPLE_FILE_03_EPH = os.path.join(dir_3, 'Spitzer_ephemeris_01.dat')  # UTC
SAMPLE_FILE_03_REF = os.path.join(dir_2, 'ob140939_Spitzer_ref_v1.dat')  # HJD'

# All in HJD':
SAMPLE_ANNUAL_PARALLAX_FILE_01 = os.path.join(dir_2, 'parallax_test_1.dat')
SAMPLE_ANNUAL_PARALLAX_FILE_02 = os.path.join(dir_2, 'parallax_test_2.dat')
SAMPLE_ANNUAL_PARALLAX_FILE_03 = os.path.join(dir_2, 'parallax_test_3.dat')
SAMPLE_ANNUAL_PARALLAX_FILE_04 = os.path.join(dir_2, 'parallax_test_4.dat')
SAMPLE_ANNUAL_PARALLAX_FILE_05 = os.path.join(dir_2, 'parallax_test_5.dat')


class _ParallaxFile(object):
    """
    Private class to allow easy access to the contents of the parallax
    test files.
    """

    def __init__(self, filename):
        """
        Open the file and store parameters.
        """
        self.filename = filename
        self.data = np.genfromtxt(
            self.filename, dtype=None,
            names=['Time', 'Magnification', 'PLflux', 'u', 'qn', 'qe'])

        (self.ulens_params, self.event_params) = self.get_file_params()

    def get_file_params(self):
        """Read in the model parameters used to create the file"""
        with open(self.filename) as data_file:
            lines = data_file.readlines()
        ulens_params = lines[3].split()
        event_params = lines[4].split()
        return (ulens_params, event_params)

    @property
    def parameters(self):
        """Model parameters"""
        model_parameters = mm.ModelParameters({
            't_0': float(self.ulens_params[1])+2450000.,
            'u_0': float(self.ulens_params[3]),
            't_E': float(self.ulens_params[4]),
            'pi_E_N': float(self.ulens_params[5]),
            'pi_E_E': float(self.ulens_params[6]),
            't_0_par': self.t_0_par})
        return model_parameters

    @property
    def coords(self):
        """Coordinates of event"""
        coords = SkyCoord(
            self.event_params[1]+' '+self.event_params[2],
            unit=(u.deg, u.deg))
        return coords

    @property
    def t_0_par(self):
        """Parallax reference time"""
        return float(self.ulens_params[2])+2450000.

    def setup_model(self):
        """Return a model using the parameters of this file"""
        model = mm.Model(parameters=self.parameters,
                         coords=self.coords)
        return model

    def setup_trajectory(self):
        """Return a trajectory using the parameters of this file"""
        trajectory = mm.Trajectory(
            self.data['Time']+2450000., parameters=self.parameters,
            parallax={'earth_orbital': True}, coords=self.coords)
        return trajectory


def test_annual_parallax_calculation():
    """
    This is a high-level unit test for parallax. The "true" values were
    calculated from the sfit routine assuming fs=1.0, fb=0.0.
    """
    t_0 = 2457479.5  # April 1 2016, a time when parallax is large
    times = np.array([t_0-1., t_0, t_0+1., t_0+1.])
    true_no_par = np.array([7.12399067, 10.0374609, 7.12399067, 7.12399067])
    true_with_par = np.array([7.12376832, 10.0386009, 7.13323363, 7.13323363])

    model_with_par = mm.Model(
        {'t_0': t_0, 'u_0': 0.1, 't_E': 10., 'pi_E': (0.3, 0.5)},
        coords='17:57:05 -30:22:59')
    model_with_par.parallax(satellite=False, earth_orbital=True,
                            topocentric=False)
    # ones = np.ones(len(times))
    # data = mm.MulensData(data_list=[times, ones, ones])
    # model_with_par.set_datasets([data])

    model_with_par.parameters.t_0_par = 2457479.

    model_no_par = mm.Model(
        {'t_0': t_0, 'u_0': 0.1, 't_E': 10.},
        coords='17:57:05 -30:22:59')
    # model_no_par.set_datasets([data])
    model_no_par.parallax(
        satellite=False, earth_orbital=False, topocentric=False)

    # Old architectures
    # np.testing.assert_almost_equal(
    #     model_no_par.data_magnification, true_no_par)
    # np.testing.assert_almost_equal(
    #     model_with_par.data_magnification, true_with_par, decimal=4)
    # New architecture
    np.testing.assert_almost_equal(
        model_no_par.get_magnification(times), true_no_par)
    np.testing.assert_almost_equal(
        model_with_par.get_magnification(times), true_with_par, decimal=4)


class test(unittest.TestCase):
    def test_wrong_settings(self):
        """
        make sure that if pi_E is defined, then at least
        one component of parallax is True
        """
        model_no_par = mm.Model(
            {'t_0': 2457479.5, 'u_0': 0.1, 't_E': 10., 'pi_E': (0.3, 0.5)},
            coords='17:57:05 -30:22:59')
        model_no_par.parallax(
            satellite=False, earth_orbital=False, topocentric=False)
        with self.assertRaises(ValueError):
            model_no_par.get_magnification(2457500.)


def do_get_delta_annual_test(filename):
    """run a test on private method Trajectory._get_delta_annual()"""
    parallax_file = _ParallaxFile(filename)
    trajectory = parallax_file.setup_trajectory()

    result = trajectory._get_delta_annual()

    np.testing.assert_almost_equal(result['N'], parallax_file.data['qn'],
                                   decimal=4)
    np.testing.assert_almost_equal(result['E'], parallax_file.data['qe'],
                                   decimal=4)


def test_get_delta_annual_1():
    """test private method Trajectory._get_delta_annual()"""
    do_get_delta_annual_test(SAMPLE_ANNUAL_PARALLAX_FILE_01)


def test_get_delta_annual_2():
    """test private method Trajectory._get_delta_annual()"""
    do_get_delta_annual_test(SAMPLE_ANNUAL_PARALLAX_FILE_02)


def test_get_delta_annual_3():
    """test private method Trajectory._get_delta_annual()"""
    do_get_delta_annual_test(SAMPLE_ANNUAL_PARALLAX_FILE_03)


def test_get_delta_annual_4():
    """test private method Trajectory._get_delta_annual()"""
    do_get_delta_annual_test(SAMPLE_ANNUAL_PARALLAX_FILE_04)


def test_get_delta_annual_5():
    """test private method Trajectory._get_delta_annual()"""
    do_get_delta_annual_test(SAMPLE_ANNUAL_PARALLAX_FILE_05)


def do_annual_parallax_test(filename):
    """testing functions called by a few unit tests"""
    with open(filename) as data_file:
        lines = data_file.readlines()
    ulens_params = lines[3].split()
    event_params = lines[4].split()
    data = np.loadtxt(filename, dtype=None)

    model = mm.Model({
        't_0': float(ulens_params[1])+2450000.,
        'u_0': float(ulens_params[3]),
        't_E': float(ulens_params[4]),
        'pi_E_N': float(ulens_params[5]),
        'pi_E_E': float(ulens_params[6])},
        coords=SkyCoord(
            event_params[1]+' '+event_params[2], unit=(u.deg, u.deg)))
    model.parameters.t_0_par = float(ulens_params[2])+2450000.

    time = data[:, 0] + 2450000.
    model.parallax(satellite=False, earth_orbital=True, topocentric=False)

    return np.testing.assert_almost_equal(
        model.get_magnification(time) / data[:, 1], 1.0, decimal=4)


def test_annual_parallax_calculation_2():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_01)


def test_annual_parallax_calculation_3():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_02)


def test_annual_parallax_calculation_4():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_03)


def test_annual_parallax_calculation_5():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_04)


def test_annual_parallax_calculation_6():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_05)


def test_satellite_and_annual_parallax_calculation():
    """test parallax calculation with Spitzer data"""
    model_parameters = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
                        'pi_E_N': -0.248, 'pi_E_E': 0.234,
                        't_0_par': 2456836.2}
    coords = "17:47:12.25 -21:22:58.2"

    model_with_par = mm.Model(model_parameters, coords=coords)
    model_with_par.parallax(satellite=True, earth_orbital=True,
                            topocentric=False)

    data_OGLE = mm.MulensData(file_name=SAMPLE_FILE_02)
    data_Spitzer = mm.MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH)
    # model_with_par.set_datasets([data_OGLE, data_Spitzer])

    model_spitzer = mm.Model(
        model_parameters, coords=coords, ephemerides_file=SAMPLE_FILE_03_EPH)

    ref_OGLE = np.loadtxt(SAMPLE_FILE_02_REF, unpack=True, usecols=[5])
    ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True, usecols=[5])

    # np.testing.assert_almost_equal(model_with_par.data_magnification[0],
    #                                ref_OGLE, decimal=2)
    # ratio = model_with_par.data_magnification[1] / ref_Spitzer
    # np.testing.assert_almost_equal(ratio, [1.]*len(ratio), decimal=4)
    np.testing.assert_almost_equal(
        model_with_par.get_magnification(data_OGLE.time), ref_OGLE,
        decimal=2)
    ratio = model_spitzer.get_magnification(data_Spitzer.time) / ref_Spitzer
    np.testing.assert_almost_equal(ratio, [1.]*len(ratio), decimal=4)


def test_satellite_parallax_magnification():
    """
    On a given date, the magnification should be different from the
    ground and from Spitzer. Use OB140939 as a test case and t0 as the
    test time.
    """
    t_0 = 2456836.22
    u_0 = 0.922
    t_E = 22.87
    pi_E_N = -0.248
    pi_E_E = 0.234

    ground_model = mm.Model(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E': [pi_E_N, pi_E_E]},
        coords='17:47:12.25 -21:22:58.2')
    space_model = mm.Model(
        {'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E': [pi_E_N, pi_E_E]},
        ra='17:47:12.25', dec='-21:22:58.2',
        ephemerides_file=SAMPLE_FILE_03_EPH)

    ground_mag = ground_model.get_magnification(t_0)
    delta = ground_mag - space_model.get_magnification(t_0)
    assert np.abs(delta) > 0.01


def test_horizons_3d():
    """
    Test if Horizons properly reads file with 3D coordinates.
    We use the satellite that has the same position as Earth.
    """
    file_in = os.path.join(dir_2, "earth_position_1.dat")
    file_out = os.path.join(dir_2, "earth_position_2.dat")

    satellite = mm.SatelliteSkyCoord(file_in)
    (times, vec_x, vec_y, vec_z) = np.loadtxt(file_out, unpack=True)
    output = satellite.get_satellite_coords(times).cartesian

    np.testing.assert_almost_equal(vec_x, output.x)
    np.testing.assert_almost_equal(vec_y, output.y)
    np.testing.assert_almost_equal(vec_z, output.z)
