from os.path import join
import numpy as np

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
