import os
import unittest

import numpy as np

import MulensModel as mm
from MulensModel.utils import PlotUtils


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


class TestGetYValueYErr(unittest.TestCase):
    """
    Tests for PlotUtils.get_y_value_y_err — addresses the
    `utils.py get_y_value_y_err()` checklist item in issue #65.
    """
    def setUp(self):
        self.flux = np.array([1.0, 10.0, 100.0])
        self.flux_err = np.array([0.01, 0.1, 1.0])

    def test_flux_format_returns_inputs_unchanged(self):
        (value, uncertainty) = PlotUtils.get_y_value_y_err(
            'flux', self.flux, self.flux_err)
        np.testing.assert_array_equal(value, self.flux)
        np.testing.assert_array_equal(uncertainty, self.flux_err)

    def test_mag_format_returns_hardcoded_magnitudes(self):
        """
        Hard-coded reference avoids defining the expected result through
        the same conversion the SUT delegates to.
        """
        (value, uncertainty) = PlotUtils.get_y_value_y_err(
            'mag', self.flux, self.flux_err)
        # MAG_ZEROPOINT = 22.0, mag = 22 - 2.5 * log10(flux)
        np.testing.assert_array_almost_equal(value, [22.0, 19.5, 17.0])
        # 2.5/ln(10) * err/flux; flux/err is constant 100 for this input
        np.testing.assert_array_almost_equal(
            uncertainty, [0.010857362, 0.010857362, 0.010857362])

    def test_unrecognized_format_raises_value_error(self):
        with self.assertRaisesRegex(
                ValueError, "Unrecognized photometry format"):
            PlotUtils.get_y_value_y_err(
                'magnitude', self.flux, self.flux_err)

    def test_uppercase_mag_is_unrecognized(self):
        """Case-sensitivity: only lowercase 'mag' / 'flux' are accepted."""
        with self.assertRaisesRegex(
                ValueError, "Unrecognized photometry format"):
            PlotUtils.get_y_value_y_err('MAG', self.flux, self.flux_err)

    def test_empty_arrays_flux_format(self):
        (value, uncertainty) = PlotUtils.get_y_value_y_err(
            'flux', np.array([]), np.array([]))
        self.assertEqual(value.size, 0)
        self.assertEqual(uncertainty.size, 0)


class TestFindSubtractXlabel(unittest.TestCase):
    """
    Tests for PlotUtils.find_subtract_xlabel — addresses the
    `utils.py find_subtract_xlabel()` checklist item in issue #65.
    """
    def test_default_returns_time(self):
        self.assertEqual(PlotUtils.find_subtract_xlabel(), 'Time')

    def test_subtract_2450000_returns_offset_label(self):
        self.assertEqual(
            PlotUtils.find_subtract_xlabel(subtract_2450000=True),
            'Time - 2450000')

    def test_subtract_2460000_returns_offset_label(self):
        self.assertEqual(
            PlotUtils.find_subtract_xlabel(subtract_2460000=True),
            'Time - 2460000')

    def test_both_subtract_flags_raise_value_error(self):
        with self.assertRaisesRegex(ValueError, "cannot be both True"):
            PlotUtils.find_subtract_xlabel(
                subtract_2450000=True, subtract_2460000=True)
