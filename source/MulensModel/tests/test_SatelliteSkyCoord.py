import unittest
from os.path import join

import numpy as np
from astropy.coordinates import SkyCoord

import MulensModel as mm


def test_boundary():
    """
    Test if calculation of position works correctly
    for epochs at the edge of range provided.
    """
    ephemeris_file = join(
        mm.DATA_PATH, 'ephemeris_files', "Spitzer_ephemeris_01.dat")

    satellite = mm.SatelliteSkyCoord(ephemeris_file)

    result_1 = satellite.get_satellite_coords([2456445.0])
    ra_1 = 15 * (8 + 26 / 60. + 37.19 / 3600.)
    dec_1 = 18 + 30 / 60. + 37.4 / 3600.
    np.testing.assert_almost_equal(result_1.ra.value, ra_1, decimal=3)
    np.testing.assert_almost_equal(result_1.dec.value, dec_1, decimal=3)

    result_2 = satellite.get_satellite_coords([2457328.0])
    ra_2 = 15 * (17 + 40 / 60. + 4.98 / 3600.)
    dec_2 = -23 - 26 / 60. - 38.2 / 3600.
    np.testing.assert_almost_equal(result_2.ra.value, ra_2, decimal=3)
    np.testing.assert_almost_equal(result_2.dec.value, dec_2, decimal=3)


class TestCheckTimes(unittest.TestCase):
    """
    Tests for SatelliteSkyCoord._check_times — addresses the
    `satelliteskycoord.py _check_times() - raise ValueError`
    checklist item in issue #65. The ephemeris file used here covers
    ~2456445 to ~2457328 (Spitzer 2013-2015).
    """
    def setUp(self):
        ephemeris_file = join(
            mm.DATA_PATH, 'ephemeris_files', 'Spitzer_ephemeris_01.dat')
        self.satellite = mm.SatelliteSkyCoord(ephemeris_file)

    def test_times_above_ephemeris_range_raise_value_error(self):
        with self.assertRaisesRegex(
                ValueError, "Satellite ephemeris doesn't cover"):
            self.satellite.get_satellite_coords([2457500.0])

    def test_times_below_ephemeris_range_raise_value_error(self):
        with self.assertRaisesRegex(
                ValueError, "Satellite ephemeris doesn't cover"):
            self.satellite.get_satellite_coords([2456000.0])

    def test_times_within_range_returns_skycoord(self):
        times = [2456800.0, 2456900.0]
        result = self.satellite.get_satellite_coords(times)
        self.assertIsInstance(result, SkyCoord)
        self.assertEqual(len(result), len(times))
