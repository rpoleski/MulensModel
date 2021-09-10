import numpy as np
import os

import MulensModel as mm


SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "photometry_files", "OB08092", "phot_ob08092_O4.dat")


def test_mag_zeropoint():
    """
    Check if changing zero-point makes global difference.
    """
    assert mm.Utils.get_flux_from_mag(22) == 1.
    original = mm.utils.MAG_ZEROPOINT
    mm.utils.MAG_ZEROPOINT = 24.5
    assert mm.Utils.get_flux_from_mag(22) == 10.
    mm.utils.MAG_ZEROPOINT = original
    assert mm.Utils.get_flux_from_mag(22) == 1.


def test_complex_fsum_1():
    z = [(0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1-1e+99j), (0.1+0.1j),
         (0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1+1e+99j), (0.1+0.1j)]
    assert mm.Utils.complex_fsum(z) == (1 + 0.8j)


def do_mag2flux_conversions_test(mag, mag_err):
    (flux, flux_err) = mm.Utils.get_flux_and_err_from_mag(mag, mag_err)
    (new_mag, new_mag_err) = mm.Utils.get_mag_and_err_from_flux(flux, flux_err)
    np.testing.assert_almost_equal(new_mag, mag, decimal=6)
    np.testing.assert_almost_equal(new_mag_err, mag_err, decimal=6)


def test_mag_and_flux_conversions_1():
    mag = 22.
    mag_err = 0.01

    (flux, flux_err) = mm.Utils.get_flux_and_err_from_mag(mag, mag_err)

    assert flux == 1.
    do_mag2flux_conversions_test(mag, mag_err)


def test_mag_and_flux_conversions_2():
    for mag in np.arange(22., 15., -1.):
        for mag_err in [0.01, 0.001, 0.1]:
            do_mag2flux_conversions_test(mag, mag_err)


def test_n_caustics():
    """test calculation of number of caustics for binary lens"""
    q = 0.123
    s_1 = 0.7605174634
    s_2 = 0.7605174635
    s_3 = 1.7289467512
    s_4 = 1.7289467513

    assert mm.Utils.get_n_caustics(s=s_1, q=q) == 3
    assert mm.Utils.get_n_caustics(s=s_2, q=q) == 1
    assert mm.Utils.get_n_caustics(s=s_3, q=q) == 1
    assert mm.Utils.get_n_caustics(s=s_4, q=q) == 2
