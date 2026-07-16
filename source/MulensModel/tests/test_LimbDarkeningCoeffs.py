import unittest
from numpy.testing import assert_almost_equal

import MulensModel as mm


class TestLimbDarkeningCoeffs(unittest.TestCase):
    """
    Unit tests for LimbDarkeningCoeffs — addresses the
    `limbdarkeningcoeffs.py` checklist item in issue #65.
    """
    def setUp(self):
        self.coeffs = mm.LimbDarkeningCoeffs()

    def test_set_and_get_gamma(self):
        self.coeffs.set_limb_coeff_gamma('I', 0.45)
        self.assertEqual(self.coeffs.get_limb_coeff_gamma('I'), 0.45)

    def test_set_and_get_u_round_trip(self):
        self.coeffs.set_limb_coeff_u('I', 0.5)
        assert_almost_equal(self.coeffs.get_limb_coeff_u('I'), 0.5)

    def test_set_u_stored_as_gamma(self):
        """set_limb_coeff_u converts via Utils.u_to_gamma before storing."""
        self.coeffs.set_limb_coeff_u('I', 0.6)
        assert_almost_equal(
            self.coeffs.get_limb_coeff_gamma('I'),
            mm.Utils.u_to_gamma(0.6))

    def test_get_gamma_unknown_bandpass_raises(self):
        with self.assertRaises(KeyError):
            self.coeffs.get_limb_coeff_gamma('Z')

    def test_get_u_unknown_bandpass_raises(self):
        with self.assertRaises(KeyError):
            self.coeffs.get_limb_coeff_u('Z')

    def test_repr_contains_set_value(self):
        self.coeffs.set_limb_coeff_gamma('I', 0.45)
        text = repr(self.coeffs)
        self.assertIn('I', text)
        self.assertIn('0.45', text)


class TestWeightedLimbCoeffGamma(unittest.TestCase):
    """
    Covers get_weighted_limb_coeff_gamma — the previously untested function
    flagged in issue #65. The pre-fix implementation omitted the per-band
    weight in the numerator (`gamma_sum += gamma` instead of
    `gamma_sum += weight * gamma`), so unequal weights silently produced
    the wrong value.
    """
    def setUp(self):
        self.coeffs = mm.LimbDarkeningCoeffs()
        self.coeffs.set_limb_coeff_gamma('I', 0.45)
        self.coeffs.set_limb_coeff_gamma('V', 0.65)

    def test_single_band(self):
        assert_almost_equal(
            self.coeffs.get_weighted_limb_coeff_gamma({'I': 1.0}), 0.45)

    def test_equal_weights(self):
        result = self.coeffs.get_weighted_limb_coeff_gamma(
            {'I': 1., 'V': 1.})
        assert_almost_equal(result, (0.45 + 0.65) / 2.)

    def test_unequal_weights_proper_weighted_mean(self):
        """
        With weights I=1.5, V=1.0 the proper weighted mean is
        (1.5*0.45 + 1.0*0.65) / (1.5 + 1.0) = 1.325 / 2.5 = 0.53.
        """
        result = self.coeffs.get_weighted_limb_coeff_gamma({'I': 1.5, 'V': 1.0})
        assert_almost_equal(result, 0.53)

    def test_non_dict_raises_type_error(self):
        with self.assertRaises(TypeError):
            self.coeffs.get_weighted_limb_coeff_gamma([('I', 1.0)])

    def test_unset_band_raises_key_error(self):
        with self.assertRaises(KeyError):
            self.coeffs.get_weighted_limb_coeff_gamma({'Z': 1.0})

    def test_empty_weights_raises_zero_division(self):
        """Document current behavior: empty dict makes weight_sum == 0."""
        with self.assertRaises(ZeroDivisionError):
            self.coeffs.get_weighted_limb_coeff_gamma({})

    def test_zero_total_weight_raises_zero_division(self):
        """Document current behavior: weights summing to 0 hit 0-division."""
        with self.assertRaises(ZeroDivisionError):
            self.coeffs.get_weighted_limb_coeff_gamma({'I': 1.0, 'V': -1.0})
