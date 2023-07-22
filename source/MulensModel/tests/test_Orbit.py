import numpy as np
from numpy.testing import assert_almost_equal
import unittest

from MulensModel.orbits.orbit import Orbit, OrbitCircular, OrbitEccentric


def test_1_circular():
    """
    circular orbit - simples calculation
    """
    orbit = OrbitCircular(period=365., semimajor_axis=1., Omega_node=0.,
                          inclination=0., ascending_node_epoch=0.)
    position = orbit.get_orbital_plane_position(time=0.)
    assert_almost_equal(position, [1., 0.])


def test_2_circular():
    """
    circular orbit at half period position
    """
    orbit = OrbitCircular(200., 1.234, 0., 0., 0.)
    position = orbit.get_orbital_plane_position(time=100.)
    assert_almost_equal(position, [-1.234, 0.])


def test_3_circular():
    """
    circular orbit, Omega = 180, and half period position
    """
    orbit = OrbitCircular(200., 1.234, 180., 0., 0.)
    position = orbit.get_reference_plane_position(time=100.)
    assert_almost_equal(position, [1.234, 0.])


def test_4_circular():
    """
    circular orbit and non-zero inclination and periapsis epoch
    """
    orbit = OrbitCircular(200., 2.345, 0., 60., 123.)
    position = orbit.get_reference_plane_position(time=123.+200/4.)
    assert_almost_equal(position, [0, 2.345/2.])


def test_5_time_vector_circular():
    """
    circular orbit calculation for multiple periods
    """
    n_repeat = 10
    orbit = OrbitCircular(200., 2.345, 0., 60., 123.)
    times = (np.arange(n_repeat) + 0.25 - 5.) * 200. + 123.
    positions = orbit.get_reference_plane_position(times)
    expected = np.repeat([[0.], [2.345/2.]], n_repeat, axis=1)
    assert_almost_equal(positions, expected)


def test_6_eccentric():
    """
    Eccentric orbit at 1 period
    """
    orbit = OrbitEccentric(
        period=400., semimajor_axis=5., Omega_node=0., inclination=0.,
        eccentricity=0.6, omega_periapsis=0., periapsis_epoch=100.)
    position = orbit.get_orbital_plane_position(300.)
    assert_almost_equal(position, [-8., 0.])


def test_7_time_vector_eccentric():
    """
    Very eccentric orbit at 1 period and half period
    """
    orbit = OrbitEccentric(400., 10., 0., 0., 0.999, 0., -100.)
    times = np.array([-500, -300])
    positions = orbit.get_orbital_plane_position(times)
    expected = np.array([[.01, -19.99], [0., 0.]])
    assert_almost_equal(positions, expected)


class TestForWrongValues(unittest.TestCase):
    def test_8_negative_period(self):
        """
        Check negative period
        """
        with self.assertRaises(ValueError):
            OrbitCircular(-365., 1., 0., 0., 123.)

    def test_9_negative_semimajor_axis(self):
        """
        Test for negative semi-major axis
        """
        with self.assertRaises(ValueError):
            OrbitCircular(365., -1., 0., 0., 123.)


def test_10_true_anomaly_large_eccentricity():
    """
    a few epochs and eccentric orbit
    """
    orbit = OrbitEccentric(400., 100., 0., 0., 0.9, 0., 0.)
    times = np.array([100., 500., 300., 150., 5.])
    true_anomalies = orbit.get_true_anomaly_deg(times)
    value = 167.70030551721663  # Value expected for input of 100 d.
    expected = np.array([
        value, value, 360.-value, 174.41306950496173, 101.28599627247006])
    assert_almost_equal(true_anomalies, expected)


def test_11_true_anomaly_huge_eccentricity():
    """
    calculate true anomaly for extremely eccentric orbit
    """
    orbit = OrbitEccentric(400., 100., 0., 0., 0.999, 0., 0.)
    true_anomaly = orbit.get_true_anomaly_deg(5.)
    assert_almost_equal(true_anomaly, 173.80472546078212, 3)


def test_12_Orbit_class_circular():
    """
    Orbit class and simplest calculation for circular orbit
    """
    orbit = Orbit(period=365., semimajor_axis=1., Omega_node=0.,
                  inclination=0., ascending_node_epoch=0.)
    position = orbit.get_orbital_plane_position(time=0.)
    assert_almost_equal(position, [1., 0.])


def test_13_Orbit_class_eccentric():
    """
    Orbit class and eccentric orbit calculation
    """
    orbit = Orbit(
        period=400., semimajor_axis=100., Omega_node=0., inclination=0.,
        eccentricity=0.9, omega_periapsis=0., periapsis_epoch=0.)
    true_anomaly = orbit.get_true_anomaly_deg(5.)
    assert_almost_equal(true_anomaly, 101.28599627247006)


class Test_Orbit_fail(unittest.TestCase):
    def test_14_Orbit_not_enough_params(self):
        """
        Orbit class with wrong input
        """
        with self.assertRaises(RuntimeError):
            Orbit(period=10., semimajor_axis=10., Omega_node=0.,
                  inclination=0.)


def test_15_OrbitCircular_based_on_argument_of_latitude():
    """
    circular orbit and non-zero argument_of_latitude_reference
    """
    orbit = OrbitCircular(
        period=365, semimajor_axis=1.5, Omega_node=12.34, inclination=-60.,
        argument_of_latitude_reference=90, epoch_reference=2456789.01234)
    position = orbit.get_orbital_plane_position(time=2456789.01234)
    assert_almost_equal(position, [0., 1.5])


class Test_OrbitCircular_fail(unittest.TestCase):
    def test_16_added_epoch(self):
        """
        missing eccentricity
        """
        with self.assertRaises(RuntimeError):
            OrbitCircular(
                period=123., semimajor_axis=5., Omega_node=90., inclination=0.,
                ascending_node_epoch=2450000., epoch_reference=2450000.)

    def test_17_added_u(self):
        """
        too many arguments for a reference epoch
        """
        with self.assertRaises(RuntimeError):
            OrbitCircular(
                period=123., semimajor_axis=5., Omega_node=90., inclination=0.,
                ascending_node_epoch=2450000.,
                argument_of_latitude_reference=90.)


def test_18_OrbitEccentric_based_on_argument_of_latitude():
    """
    Eccentric orbit calculation starting from
    argument_of_latitude_reference
    """
    orbit = OrbitEccentric(
        period=365, semimajor_axis=1.5, Omega_node=12.34, inclination=-60.,
        eccentricity=0.5, omega_periapsis=0.,
        argument_of_latitude_reference=0, epoch_reference=2456789.01234)
    position = orbit.get_orbital_plane_position(2456789.01234)
    assert_almost_equal(position, [0.75, 0.])


def test_19_OrbitEccentric_based_on_argument_of_latitude():
    """
    Eccentric orbit calculation starting from non-zero value of
    argument_of_latitude_reference
    """
    orbit = OrbitEccentric(
        period=360., semimajor_axis=1.5, Omega_node=12.34, inclination=-60.,
        eccentricity=0.5,
        omega_periapsis=0., epoch_reference=2456789.01234,
        argument_of_latitude_reference=118.81500092699673)
    # An independend code says the above value is nu for e=0.5 and t=t_0+P/6.
    position = orbit.get_orbital_plane_position(2456789.01234-60.)
    assert_almost_equal(position, [0.75, 0.])


def test_20_OrbitEccentric_based_on_argument_of_latitude():
    """
    Eccentric orbit with all argument non-zero and
    argument_of_latitude_reference in __init__()
    """
    orbit = OrbitEccentric(
        period=360., semimajor_axis=1.5, Omega_node=12.34, inclination=-60.,
        eccentricity=0.5,
        omega_periapsis=100.,
        argument_of_latitude_reference=118.81500092699673+100,
        epoch_reference=2456789.01234+60)
    position = orbit.get_orbital_plane_position(2456789.01234-180.)
    assert_almost_equal(position, [-2.25, 0.])
