import warnings
import numpy as np
from math import sqrt

from MulensModel.binarylens import BinaryLens
from MulensModel.binarylensimports import (_vbbl_wrapped,
                                           _vbbl_binary_mag_0, _vbbl_SG12_9)
import MulensModel as mm


class BinaryLensWithShear(BinaryLens):
    """
    The binary lens equation - its solutions, images, parities,
    magnifications, etc.

    The binary lens equation is a 5th order complex polynomial
    or a 9th order complex polynomial if including external shear.

    Attributes :
        mass_1: *float*
            mass of the primary (left-hand object) as a fraction of
            the total mass.

        mass_2: *float*
            mass of the secondary (right-hand object) as a fraction of the
            total mass.

        separation: *float*
            separation between the two bodies as a fraction of the Einstein
            ring.

        convergence_K: *float*
            External mass sheet convergence.

        shear_G: *complex*
            External mass sheat shear.

    Note: mass_1 and mass_2 may be defined as a fraction of some other
    mass than the total mass. This is possible but not recommended -
    make sure you know what you're doing before you start using this
    possibility.

    If you're using this class, then please cite
    Peirson et al. (2022; ApJ 927, 24).
    """

    def __init__(self, mass_1=None, mass_2=None, separation=None,
                 convergence_K=None, shear_G=None):
        BinaryLens.__init__(self, mass_1, mass_2, separation)
        self.convergence_K = convergence_K
        self.shear_G = shear_G

    def _get_polynomial_WM95(self, source_x, source_y):
        """
        calculate coefficients of the polynomial in geometric center frame
        """
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff

        pos_z1 = self._position_z1

        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2
        z1_pow5 = z1_pow2 * z1_pow2 * pos_z1
        z1_pow6 = z1_pow2 * z1_pow2 * z1_pow2

        # Convergence added here
        convergence_K = self.convergence_K
        K_pow2 = convergence_K * convergence_K
        K_pow3 = convergence_K * K_pow2

        zeta = self._zeta
        zeta_conj = zeta.conjugate()
        zeta_conj_pow2 = zeta_conj * zeta_conj

        # Calculate the coefficients of the 5th order complex polynomial now
        # with external convergence
        coeff_5 = mm.Utils.complex_fsum(
            [-z1_pow2 * K_pow3, 3 * z1_pow2 * K_pow2, -3 * z1_pow2 *
             convergence_K, z1_pow2, convergence_K * zeta_conj_pow2, -
             zeta_conj_pow2])
        coeff_4 = mm.Utils.complex_fsum(
            [2 * total_m * convergence_K *
             zeta_conj, -2 * total_m * zeta_conj, -
             zeta * z1_pow2 * K_pow2, 2 * zeta *
             z1_pow2 * convergence_K, -zeta *
             z1_pow2, zeta * zeta_conj_pow2, 2 *
             pos_z1 * K_pow2 * m_diff, -4 * pos_z1
             * convergence_K * m_diff, 2 * pos_z1 *
             m_diff])
        coeff_3 = mm.Utils.complex_fsum(
            [4 * total_m * zeta * zeta_conj, 2 * z1_pow4 * K_pow3, -6 * z1_pow4
             * K_pow2, 6 * z1_pow4 * convergence_K, -2 * z1_pow4, -2 * z1_pow2
             * convergence_K * zeta_conj_pow2, 2 * z1_pow2 * zeta_conj_pow2, 4
             * pos_z1 * convergence_K * m_diff * zeta_conj, -4 * pos_z1 *
             m_diff * zeta_conj])
        coeff_2 = mm.Utils.complex_fsum(
            [4 * total_m_pow2 * zeta, 4 * total_m * pos_z1 * convergence_K *
             m_diff, -4 * total_m * pos_z1 * m_diff, 2 * zeta * z1_pow4 *
             K_pow2, -4 * zeta * z1_pow4 * convergence_K, 2 * zeta * z1_pow4, -
             2 * zeta * z1_pow2 * zeta_conj_pow2, 4 * zeta * pos_z1 * m_diff *
             zeta_conj, -4 * z1_pow3 * K_pow2 * m_diff, 8 * z1_pow3 *
             convergence_K * m_diff, -4 * z1_pow3 * m_diff])
        coeff_1 = mm.Utils.complex_fsum(
            [4 * total_m_pow2 * z1_pow2 * convergence_K, -4 * total_m_pow2 *
             z1_pow2, -4 * total_m * zeta * z1_pow2 * zeta_conj, 8 * total_m *
             zeta * pos_z1 * m_diff, -z1_pow6 * K_pow3, 3 * z1_pow6 * K_pow2, -
             3 * z1_pow6 * convergence_K, z1_pow6, z1_pow4 * convergence_K *
             zeta_conj_pow2, -z1_pow4 * zeta_conj_pow2, -4 * z1_pow3 *
             convergence_K * m_diff * zeta_conj, 4 * z1_pow3 * m_diff *
             zeta_conj, 4 * z1_pow2 * convergence_K * m_diff_pow2, -4 * z1_pow2
             * m_diff_pow2])
        coeff_0 = mm.Utils.complex_fsum(
            [-2 * total_m * z1_pow4 * convergence_K * zeta_conj, 2 * total_m *
             z1_pow4 * zeta_conj, 4 * total_m * z1_pow3 * convergence_K *
             m_diff, -4 * total_m * z1_pow3 * m_diff, -zeta * z1_pow6 * K_pow2,
             2 * zeta * z1_pow6 * convergence_K, -zeta * z1_pow6, zeta *
             z1_pow4 * zeta_conj_pow2, -4 * zeta * z1_pow3 * m_diff *
             zeta_conj, 4 * zeta * z1_pow2 * m_diff_pow2, 2 * z1_pow5 * K_pow2
             * m_diff, -4 * z1_pow5 * convergence_K * m_diff, 2 * z1_pow5 *
             m_diff])

        # Return the coefficients of the polynomial
        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5]
        return np.array(coeffs_list).reshape(6)

    def _get_polynomial_planet_frame(self, source_x, source_y):
        """calculate coefficients of the polynomial in planet frame"""
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m
        total_m_pow3 = total_m * total_m_pow2

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff
        m_diff_pow3 = m_diff * m_diff_pow2

        zeta = self._zeta
        zeta_conj = zeta.conjugate()
        zeta_conj_pow2 = zeta_conj * zeta_conj
        zeta_conj_pow3 = zeta_conj * zeta_conj_pow2

        # External Convergence added here
        convergence_K = self.convergence_K
        K_pow2 = convergence_K * convergence_K
        K_pow3 = convergence_K * K_pow2

        # External Shear added here
        shear_G = self.shear_G
        Gc = np.conjugate(shear_G)

        Gc_pow2 = Gc * Gc
        Gc_pow3 = Gc * Gc_pow2

        c_sum = mm.Utils.complex_fsum

        z1 = self._position_z1
        z1_pow2 = z1 * z1
        z1_pow3 = z1_pow2 * z1
        z1_pow4 = z1_pow2 * z1_pow2

        coeff_9 = c_sum([-shear_G*Gc_pow3, Gc_pow2*K_pow2, -
                        2*Gc_pow2*convergence_K, Gc_pow2])
        coeff_8 = c_sum(
            [zeta * Gc_pow2 * convergence_K, -zeta * Gc_pow2, 3 * z1 * shear_G
             * Gc_pow3, -z1 * shear_G * Gc_pow2 * convergence_K, z1 * shear_G *
             Gc_pow2, -3 * z1 * Gc_pow2 * K_pow2, 6 * z1 * Gc_pow2 *
             convergence_K, -3 * z1 * Gc_pow2, z1 * Gc * K_pow3, -3 * z1 * Gc *
             K_pow2, 3 * z1 * Gc * convergence_K, -z1 * Gc, -3 * shear_G *
             Gc_pow2 * zeta_conj, 2 * Gc * K_pow2 * zeta_conj, -4 * Gc *
             convergence_K * zeta_conj, 2 * Gc * zeta_conj])
        coeff_7 = c_sum(
            [-6 * total_m * shear_G * Gc_pow2, 2 * total_m * Gc *
             K_pow2, -4 * total_m * Gc * convergence_K, 2 * total_m
             * Gc, -3 * zeta * z1 * Gc_pow2 * convergence_K, 3 *
             zeta * z1 * Gc_pow2, zeta * z1 * Gc * K_pow2, -2 *
             zeta * z1 * Gc * convergence_K, zeta * z1 * Gc, 2 *
             zeta * Gc * convergence_K * zeta_conj, -2 * zeta * Gc
             * zeta_conj, -3 * z1_pow2 * shear_G * Gc_pow3, 3 *
             z1_pow2 * shear_G * Gc_pow2 * convergence_K, -3 *
             z1_pow2 * shear_G * Gc_pow2, 3 * z1_pow2 * Gc_pow2 *
             K_pow2, -6 * z1_pow2 * Gc_pow2 * convergence_K, 3 *
             z1_pow2 * Gc_pow2, -3 * z1_pow2 * Gc * K_pow3, 9 *
             z1_pow2 * Gc * K_pow2, -9 * z1_pow2 * Gc *
             convergence_K, 3 * z1_pow2 * Gc, 9 * z1 * shear_G *
             Gc_pow2 * zeta_conj, -2 * z1 * shear_G * Gc *
             convergence_K * zeta_conj, 2 * z1 * shear_G * Gc *
             zeta_conj, -6 * z1 * Gc * K_pow2 * zeta_conj, 12 * z1
             * Gc * convergence_K * zeta_conj, -6 * z1 * Gc *
             zeta_conj, z1 * K_pow3 * zeta_conj, -3 * z1 * K_pow2 *
             zeta_conj, 3 * z1 * convergence_K * zeta_conj, -z1 *
             zeta_conj, -3 * shear_G * Gc * zeta_conj_pow2, K_pow2
             * zeta_conj_pow2, -2 * convergence_K * zeta_conj_pow2,
             zeta_conj_pow2])
        coeff_6 = c_sum(
            [4 * total_m * zeta * Gc * convergence_K, -4 * total_m * zeta * Gc,
             15 * total_m * z1 * shear_G * Gc_pow2, -4 * total_m * z1 * shear_G
             * Gc * convergence_K, 4 * total_m * z1 * shear_G * Gc, -4 *
             total_m * z1 * Gc * K_pow2, 8 * total_m * z1 * Gc * convergence_K,
             -4 * total_m * z1 * Gc, total_m * z1 * K_pow3, -3 * total_m * z1 *
             K_pow2, 3 * total_m * z1 * convergence_K, -total_m * z1, -12 *
             total_m * shear_G * Gc * zeta_conj, 2 * total_m * K_pow2 *
             zeta_conj, -4 * total_m * convergence_K * zeta_conj, 2 * total_m *
             zeta_conj, 3 * zeta * z1_pow2 * Gc_pow2 * convergence_K, -3 * zeta
             * z1_pow2 * Gc_pow2, -3 * zeta * z1_pow2 * Gc * K_pow2, 6 * zeta *
             z1_pow2 * Gc * convergence_K, -3 * zeta * z1_pow2 * Gc, -6 * zeta
             * z1 * Gc * convergence_K * zeta_conj, 6 * zeta * z1 * Gc *
             zeta_conj, zeta * z1 * K_pow2 * zeta_conj, -2 * zeta * z1 *
             convergence_K * zeta_conj, zeta * z1 * zeta_conj, zeta *
             convergence_K * zeta_conj_pow2, -zeta * zeta_conj_pow2, z1_pow3 *
             shear_G * Gc_pow3, -3 * z1_pow3 * shear_G * Gc_pow2 *
             convergence_K, 3 * z1_pow3 * shear_G * Gc_pow2, -z1_pow3 * Gc_pow2
             * K_pow2, 2 * z1_pow3 * Gc_pow2 * convergence_K, -z1_pow3 *
             Gc_pow2, 3 * z1_pow3 * Gc * K_pow3, -9 * z1_pow3 * Gc * K_pow2, 9
             * z1_pow3 * Gc * convergence_K, -3 * z1_pow3 * Gc, -9 * z1_pow2 *
             shear_G * Gc_pow2 * zeta_conj, 6 * z1_pow2 * shear_G * Gc *
             convergence_K * zeta_conj, -6 * z1_pow2 * shear_G * Gc *
             zeta_conj, 6 * z1_pow2 * Gc * K_pow2 * zeta_conj, -12 * z1_pow2 *
             Gc * convergence_K * zeta_conj, 6 * z1_pow2 * Gc * zeta_conj, -3 *
             z1_pow2 * K_pow3 * zeta_conj, 9 * z1_pow2 * K_pow2 * zeta_conj, -
             9 * z1_pow2 * convergence_K * zeta_conj, 3 * z1_pow2 * zeta_conj,
             3 * z1 * shear_G * Gc_pow2 * m_diff, 9 * z1 * shear_G * Gc *
             zeta_conj_pow2, -z1 * shear_G * convergence_K * zeta_conj_pow2, z1
             * shear_G * zeta_conj_pow2, -2 * z1 * Gc * K_pow2 * m_diff, 4 * z1
             * Gc * convergence_K * m_diff, -2 * z1 * Gc * m_diff, -z1 * K_pow3
             * m_diff, 3 * z1 * K_pow2 * m_diff, -3 * z1 * K_pow2 *
             zeta_conj_pow2, -3 * z1 * convergence_K * m_diff, 6 * z1 *
             convergence_K * zeta_conj_pow2, z1 * m_diff, -3 * z1 *
             zeta_conj_pow2, -shear_G * zeta_conj_pow3])
        coeff_5 = c_sum(
            [-12 * total_m_pow2 * shear_G * Gc, -10 * total_m * zeta * z1 * Gc
             * convergence_K, 10 * total_m * zeta * z1 * Gc, 2 * total_m * zeta
             * z1 * K_pow2, -4 * total_m * zeta * z1 * convergence_K, 2 *
             total_m * zeta * z1, 4 * total_m * zeta * convergence_K *
             zeta_conj, -4 * total_m * zeta * zeta_conj, -12 * total_m *
             z1_pow2 * shear_G * Gc_pow2, 10 * total_m * z1_pow2 * shear_G * Gc
             * convergence_K, -10 * total_m * z1_pow2 * shear_G * Gc, 2 *
             total_m * z1_pow2 * Gc * K_pow2, -4 * total_m * z1_pow2 * Gc *
             convergence_K, 2 * total_m * z1_pow2 * Gc, -2 * total_m * z1_pow2
             * K_pow3, 6 * total_m * z1_pow2 * K_pow2, -6 * total_m * z1_pow2 *
             convergence_K, 2 * total_m * z1_pow2, 30 * total_m * z1 * shear_G
             * Gc * zeta_conj, -4 * total_m * z1 * shear_G * convergence_K *
             zeta_conj, 4 * total_m * z1 * shear_G * zeta_conj, -4 * total_m *
             z1 * K_pow2 * zeta_conj, 8 * total_m * z1 * convergence_K *
             zeta_conj, -4 * total_m * z1 * zeta_conj, -6 * total_m * shear_G *
             zeta_conj_pow2, -zeta * z1_pow3 * Gc_pow2 * convergence_K, zeta *
             z1_pow3 * Gc_pow2, 3 * zeta * z1_pow3 * Gc * K_pow2, -6 * zeta *
             z1_pow3 * Gc * convergence_K, 3 * zeta * z1_pow3 * Gc, 6 * zeta *
             z1_pow2 * Gc * convergence_K * zeta_conj, -6 * zeta * z1_pow2 * Gc
             * zeta_conj, -3 * zeta * z1_pow2 * K_pow2 * zeta_conj, 6 * zeta *
             z1_pow2 * convergence_K * zeta_conj, -3 * zeta * z1_pow2 *
             zeta_conj, -2 * zeta * z1 * Gc * convergence_K * m_diff, 2 * zeta
             * z1 * Gc * m_diff, -3 * zeta * z1 * convergence_K *
             zeta_conj_pow2, 3 * zeta * z1 * zeta_conj_pow2, z1_pow4 * shear_G
             * Gc_pow2 * convergence_K, -z1_pow4 * shear_G * Gc_pow2, -z1_pow4
             * Gc * K_pow3, 3 * z1_pow4 * Gc * K_pow2, -3 * z1_pow4 * Gc *
             convergence_K, z1_pow4 * Gc, 3 * z1_pow3 * shear_G * Gc_pow2 *
             zeta_conj, -6 * z1_pow3 * shear_G * Gc * convergence_K *
             zeta_conj, 6 * z1_pow3 * shear_G * Gc * zeta_conj, -2 * z1_pow3 *
             Gc * K_pow2 * zeta_conj, 4 * z1_pow3 * Gc * convergence_K *
             zeta_conj, -2 * z1_pow3 * Gc * zeta_conj, 3 * z1_pow3 * K_pow3 *
             zeta_conj, -9 * z1_pow3 * K_pow2 * zeta_conj, 9 * z1_pow3 *
             convergence_K * zeta_conj, -3 * z1_pow3 * zeta_conj, -6 * z1_pow2
             * shear_G * Gc_pow2 * m_diff, 2 * z1_pow2 * shear_G * Gc *
             convergence_K * m_diff, -2 * z1_pow2 * shear_G * Gc * m_diff, -9 *
             z1_pow2 * shear_G * Gc * zeta_conj_pow2, 3 * z1_pow2 * shear_G *
             convergence_K * zeta_conj_pow2, -3 * z1_pow2 * shear_G *
             zeta_conj_pow2, 4 * z1_pow2 * Gc * K_pow2 * m_diff, -8 * z1_pow2 *
             Gc * convergence_K * m_diff, 4 * z1_pow2 * Gc * m_diff, 2 *
             z1_pow2 * K_pow3 * m_diff, -6 * z1_pow2 * K_pow2 * m_diff, 3 *
             z1_pow2 * K_pow2 * zeta_conj_pow2, 6 * z1_pow2 * convergence_K *
             m_diff, -6 * z1_pow2 * convergence_K * zeta_conj_pow2, -2 *
             z1_pow2 * m_diff, 3 * z1_pow2 * zeta_conj_pow2, 6 * z1 * shear_G *
             Gc * m_diff * zeta_conj, 3 * z1 * shear_G * zeta_conj_pow3, -2 *
             z1 * K_pow2 * m_diff * zeta_conj, 4 * z1 * convergence_K * m_diff
             * zeta_conj, -2 * z1 * m_diff * zeta_conj])
        coeff_4 = c_sum(
            [4 * total_m_pow2 * zeta * convergence_K, -4 * total_m_pow2 * zeta,
             24 * total_m_pow2 * z1 * shear_G * Gc, -4 * total_m_pow2 * z1 *
             shear_G * convergence_K, 4 * total_m_pow2 * z1 * shear_G, 2 *
             total_m_pow2 * z1 * K_pow2, -4 * total_m_pow2 * z1 *
             convergence_K, 2 * total_m_pow2 * z1, -12 * total_m_pow2 *
             shear_G * zeta_conj, 8 * total_m * zeta * z1_pow2 * Gc *
             convergence_K, -8 * total_m * zeta * z1_pow2 * Gc, -5 * total_m *
             zeta * z1_pow2 * K_pow2, 10 * total_m * zeta * z1_pow2 *
             convergence_K, -5 * total_m * zeta * z1_pow2, -10 * total_m * zeta
             * z1 * convergence_K * zeta_conj, 10 * total_m * zeta * z1 *
             zeta_conj, 3 * total_m * z1_pow3 * shear_G * Gc_pow2, -8 * total_m
             * z1_pow3 * shear_G * Gc * convergence_K, 8 * total_m * z1_pow3 *
             shear_G * Gc, total_m * z1_pow3 * K_pow3, -3 * total_m * z1_pow3 *
             K_pow2, 3 * total_m * z1_pow3 * convergence_K, -total_m * z1_pow3,
             -24 * total_m * z1_pow2 * shear_G * Gc * zeta_conj, 10 * total_m *
             z1_pow2 * shear_G * convergence_K * zeta_conj, -10 * total_m *
             z1_pow2 * shear_G * zeta_conj, 2 * total_m * z1_pow2 * K_pow2 *
             zeta_conj, -4 * total_m * z1_pow2 * convergence_K * zeta_conj, 2 *
             total_m * z1_pow2 * zeta_conj, 12 * total_m * z1 * shear_G * Gc *
             m_diff, 15 * total_m * z1 * shear_G * zeta_conj_pow2, -2 * total_m
             * z1 * K_pow2 * m_diff, 4 * total_m * z1 * convergence_K * m_diff,
             -2 * total_m * z1 * m_diff, -zeta * z1_pow4 * Gc * K_pow2, 2 *
             zeta * z1_pow4 * Gc * convergence_K, -zeta * z1_pow4 * Gc, -2 *
             zeta * z1_pow3 * Gc * convergence_K * zeta_conj, 2 * zeta *
             z1_pow3 * Gc * zeta_conj, 3 * zeta * z1_pow3 * K_pow2 * zeta_conj,
             -6 * zeta * z1_pow3 * convergence_K * zeta_conj, 3 * zeta *
             z1_pow3 * zeta_conj, 4 * zeta * z1_pow2 * Gc * convergence_K *
             m_diff, -4 * zeta * z1_pow2 * Gc * m_diff, -zeta * z1_pow2 *
             K_pow2 * m_diff, 2 * zeta * z1_pow2 * convergence_K * m_diff, 3 *
             zeta * z1_pow2 * convergence_K * zeta_conj_pow2, -zeta * z1_pow2 *
             m_diff, -3 * zeta * z1_pow2 * zeta_conj_pow2, -2 * zeta * z1 *
             convergence_K * m_diff * zeta_conj, 2 * zeta * z1 * m_diff *
             zeta_conj, 2 * z1_pow4 * shear_G * Gc * convergence_K * zeta_conj,
             -2 * z1_pow4 * shear_G * Gc * zeta_conj, -z1_pow4 * K_pow3 *
             zeta_conj, 3 * z1_pow4 * K_pow2 * zeta_conj, -3 * z1_pow4 *
             convergence_K * zeta_conj, z1_pow4 * zeta_conj, 3 * z1_pow3 *
             shear_G * Gc_pow2 * m_diff, -4 * z1_pow3 * shear_G * Gc *
             convergence_K * m_diff, 4 * z1_pow3 * shear_G * Gc * m_diff, 3 *
             z1_pow3 * shear_G * Gc * zeta_conj_pow2, -3 * z1_pow3 * shear_G *
             convergence_K * zeta_conj_pow2, 3 * z1_pow3 * shear_G *
             zeta_conj_pow2, -2 * z1_pow3 * Gc * K_pow2 * m_diff, 4 * z1_pow3 *
             Gc * convergence_K * m_diff, -2 * z1_pow3 * Gc * m_diff, -z1_pow3
             * K_pow3 * m_diff, 3 * z1_pow3 * K_pow2 * m_diff, -z1_pow3 *
             K_pow2 * zeta_conj_pow2, -3 * z1_pow3 * convergence_K * m_diff, 2
             * z1_pow3 * convergence_K * zeta_conj_pow2, z1_pow3 * m_diff, -
             z1_pow3 * zeta_conj_pow2, -12 * z1_pow2 * shear_G * Gc * m_diff *
             zeta_conj, 2 * z1_pow2 * shear_G * convergence_K * m_diff *
             zeta_conj, -2 * z1_pow2 * shear_G * m_diff * zeta_conj, -3 *
             z1_pow2 * shear_G * zeta_conj_pow3, 4 * z1_pow2 * K_pow2 * m_diff
             * zeta_conj, -8 * z1_pow2 * convergence_K * m_diff * zeta_conj, 4
             * z1_pow2 * m_diff * zeta_conj, 3 * z1 * shear_G * m_diff *
             zeta_conj_pow2])
        coeff_3 = c_sum(
            [-8 * total_m_pow3 * shear_G, -8 * total_m_pow2 * zeta * z1 *
             convergence_K, 8 * total_m_pow2 * zeta * z1, -15 * total_m_pow2 *
             z1_pow2 * shear_G * Gc, 8 * total_m_pow2 * z1_pow2 * shear_G *
             convergence_K, -8 * total_m_pow2 * z1_pow2 * shear_G, -3 *
             total_m_pow2 * z1_pow2 * K_pow2, 6 * total_m_pow2 * z1_pow2 *
             convergence_K, -3 * total_m_pow2 * z1_pow2, 24 * total_m_pow2 * z1
             * shear_G * zeta_conj, -2 * total_m * zeta * z1_pow3 * Gc *
             convergence_K, 2 * total_m * zeta * z1_pow3 * Gc, 4 * total_m *
             zeta * z1_pow3 * K_pow2, -8 * total_m * zeta * z1_pow3 *
             convergence_K, 4 * total_m * zeta * z1_pow3, 8 * total_m * zeta *
             z1_pow2 * convergence_K * zeta_conj, -8 * total_m * zeta * z1_pow2
             * zeta_conj, -4 * total_m * zeta * z1 * convergence_K * m_diff, 4
             * total_m * zeta * z1 * m_diff, 2 * total_m * z1_pow4 * shear_G *
             Gc * convergence_K, -2 * total_m * z1_pow4 * shear_G * Gc, 6 *
             total_m * z1_pow3 * shear_G * Gc * zeta_conj, -8 * total_m *
             z1_pow3 * shear_G * convergence_K * zeta_conj, 8 * total_m *
             z1_pow3 * shear_G * zeta_conj, -18 * total_m * z1_pow2 * shear_G *
             Gc * m_diff, 4 * total_m * z1_pow2 * shear_G * convergence_K *
             m_diff, -4 * total_m * z1_pow2 * shear_G * m_diff, -12 * total_m *
             z1_pow2 * shear_G * zeta_conj_pow2, 2 * total_m * z1_pow2 * K_pow2
             * m_diff, -4 * total_m * z1_pow2 * convergence_K * m_diff, 2 *
             total_m * z1_pow2 * m_diff, 12 * total_m * z1 * shear_G * m_diff *
             zeta_conj, -zeta * z1_pow4 * K_pow2 * zeta_conj, 2 * zeta *
             z1_pow4 * convergence_K * zeta_conj, -zeta * z1_pow4 * zeta_conj,
             -2 * zeta * z1_pow3 * Gc * convergence_K * m_diff, 2 * zeta *
             z1_pow3 * Gc * m_diff, 2 * zeta * z1_pow3 * K_pow2 * m_diff, -4 *
             zeta * z1_pow3 * convergence_K * m_diff, -zeta * z1_pow3 *
             convergence_K * zeta_conj_pow2, 2 * zeta * z1_pow3 * m_diff, zeta
             * z1_pow3 * zeta_conj_pow2, 4 * zeta * z1_pow2 * convergence_K *
             m_diff * zeta_conj, -4 * zeta * z1_pow2 * m_diff * zeta_conj, 2 *
             z1_pow4 * shear_G * Gc * convergence_K * m_diff, -2 * z1_pow4 *
             shear_G * Gc * m_diff, z1_pow4 * shear_G * convergence_K *
             zeta_conj_pow2, -z1_pow4 * shear_G * zeta_conj_pow2, 6 * z1_pow3 *
             shear_G * Gc * m_diff * zeta_conj, -4 * z1_pow3 * shear_G *
             convergence_K * m_diff * zeta_conj, 4 * z1_pow3 * shear_G * m_diff
             * zeta_conj, z1_pow3 * shear_G * zeta_conj_pow3, -2 * z1_pow3 *
             K_pow2 * m_diff * zeta_conj, 4 * z1_pow3 * convergence_K * m_diff
             * zeta_conj, -2 * z1_pow3 * m_diff * zeta_conj, -3 * z1_pow2 *
             shear_G * Gc * m_diff_pow2, -6 * z1_pow2 * shear_G * m_diff *
             zeta_conj_pow2, z1_pow2 * K_pow2 * m_diff_pow2, -2 * z1_pow2 *
             convergence_K * m_diff_pow2, z1_pow2 * m_diff_pow2])
        coeff_2 = c_sum(
            [12 * total_m_pow3 * z1 * shear_G, 5 * total_m_pow2 * zeta *
             z1_pow2 * convergence_K, -5 * total_m_pow2 * zeta * z1_pow2, 3 *
             total_m_pow2 * z1_pow3 * shear_G * Gc, -5 * total_m_pow2 * z1_pow3
             * shear_G * convergence_K, 5 * total_m_pow2 * z1_pow3 * shear_G,
             total_m_pow2 * z1_pow3 * K_pow2, -2 * total_m_pow2 * z1_pow3 *
             convergence_K, total_m_pow2 * z1_pow3, -15 * total_m_pow2 *
             z1_pow2 * shear_G * zeta_conj, 12 * total_m_pow2 * z1 * shear_G *
             m_diff, -total_m * zeta * z1_pow4 * K_pow2, 2 * total_m * zeta *
             z1_pow4 * convergence_K, -total_m * zeta * z1_pow4, -2 * total_m *
             zeta * z1_pow3 * convergence_K * zeta_conj, 2 * total_m * zeta *
             z1_pow3 * zeta_conj, 6 * total_m * zeta * z1_pow2 * convergence_K
             * m_diff, -6 * total_m * zeta * z1_pow2 * m_diff, 2 * total_m *
             z1_pow4 * shear_G * convergence_K * zeta_conj, -2 * total_m *
             z1_pow4 * shear_G * zeta_conj, 6 * total_m * z1_pow3 * shear_G *
             Gc * m_diff, -6 * total_m * z1_pow3 * shear_G * convergence_K *
             m_diff, 6 * total_m * z1_pow3 * shear_G * m_diff, 3 * total_m *
             z1_pow3 * shear_G * zeta_conj_pow2, -18 * total_m * z1_pow2 *
             shear_G * m_diff * zeta_conj, -zeta * z1_pow4 * K_pow2 * m_diff, 2
             * zeta * z1_pow4 * convergence_K * m_diff, -zeta * z1_pow4 *
             m_diff, -2 * zeta * z1_pow3 * convergence_K * m_diff * zeta_conj,
             2 * zeta * z1_pow3 * m_diff * zeta_conj, zeta * z1_pow2 *
             convergence_K * m_diff_pow2, -zeta * z1_pow2 * m_diff_pow2, 2 *
             z1_pow4 * shear_G * convergence_K * m_diff * zeta_conj, -2 *
             z1_pow4 * shear_G * m_diff * zeta_conj, 3 * z1_pow3 * shear_G * Gc
             * m_diff_pow2, -z1_pow3 * shear_G * convergence_K * m_diff_pow2,
             z1_pow3 * shear_G * m_diff_pow2, 3 * z1_pow3 * shear_G * m_diff *
             zeta_conj_pow2, -z1_pow3 * K_pow2 * m_diff_pow2, 2 * z1_pow3 *
             convergence_K * m_diff_pow2, -z1_pow3 * m_diff_pow2, -3 * z1_pow2
             * shear_G * m_diff_pow2 * zeta_conj])
        coeff_1 = c_sum(
            [-6 * total_m_pow3 * z1_pow2 * shear_G, -total_m_pow2 * zeta *
             z1_pow3 * convergence_K, total_m_pow2 * zeta * z1_pow3,
             total_m_pow2 * z1_pow4 * shear_G * convergence_K, -total_m_pow2 *
             z1_pow4 * shear_G, 3 * total_m_pow2 * z1_pow3 * shear_G *
             zeta_conj, -12 * total_m_pow2 * z1_pow2 * shear_G * m_diff, -2 *
             total_m * zeta * z1_pow3 * convergence_K * m_diff, 2 * total_m *
             zeta * z1_pow3 * m_diff, 2 * total_m * z1_pow4 * shear_G *
             convergence_K * m_diff, -2 * total_m * z1_pow4 * shear_G * m_diff,
             6 * total_m * z1_pow3 * shear_G * m_diff * zeta_conj, -6 * total_m
             * z1_pow2 * shear_G * m_diff_pow2, -zeta * z1_pow3 * convergence_K
             * m_diff_pow2, zeta * z1_pow3 * m_diff_pow2, z1_pow4 * shear_G *
             convergence_K * m_diff_pow2, -z1_pow4 * shear_G * m_diff_pow2, 3 *
             z1_pow3 * shear_G * m_diff_pow2 * zeta_conj])
        coeff_0 = c_sum(
            [total_m_pow3 * z1_pow3 * shear_G, 3 * total_m_pow2 * z1_pow3 *
             shear_G * m_diff, 3 * total_m * z1_pow3 * shear_G * m_diff_pow2,
             z1_pow3 * shear_G * m_diff_pow3])

        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5,
                       coeff_6, coeff_7, coeff_8, coeff_9]
        return np.array(coeffs_list).reshape(10)

    def _get_polynomial_roots(self, source_x, source_y):
        """roots of the polynomial"""
        polynomial_input = [
            self.mass_1, self.mass_2, self.separation, self.convergence_K,
            self.shear_G, source_x, source_y]

        if polynomial_input == self._last_polynomial_input:
            return self._polynomial_roots

        polynomial = self._get_polynomial(
            source_x=source_x, source_y=source_y)

        np_polyroots = np.polynomial.polynomial.polyroots
        if self._solver == 'numpy':
            self._polynomial_roots = np_polyroots(polynomial)
        elif self._solver == 'Skowron_and_Gould_12':
            args = polynomial.real.tolist() + polynomial.imag.tolist()
            try:
                out = _vbbl_SG12_9(*args)
            except ValueError as err:
                err2 = "\n\nSwitching from Skowron & Gould 2012 to numpy"
                warnings.warn(str(err) + err2, UserWarning)
                self._solver = 'numpy'
                self._polynomial_roots = np_polyroots(polynomial)
            else:
                roots = [
                    out[i] + out[i+9] * 1.j
                    for i in range(9) if (
                        (abs(out[i]) > 1e-10 or abs(out[i+9]) > 1e-10)
                        and (abs(out[i] - self._position_z1) > 1e-10
                             or abs(out[i+9]) > 1e-10))]
                self._polynomial_roots = np.array(roots)
        else:
            raise ValueError('Unknown solver: {:}'.format(self._solver))
        self._last_polynomial_input = polynomial_input

        return self._polynomial_roots

    def _verify_polynomial_roots(
            self, source_x, source_y, return_distances=False):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots(
            source_x=source_x, source_y=source_y)

        roots_conj = np.conjugate(roots)
        component2 = self.mass_1 / (roots_conj - self._position_z1)
        component3 = self.mass_2 / (roots_conj - self._position_z2)
        solutions = (self._zeta + self.shear_G * roots_conj +
                     component2 + component3) / (1 - self.convergence_K)

        # This backs-up the lens equation.

        out = []
        distances = []
        for (i, root) in enumerate(roots):
            distances_from_root = abs((solutions-root)**2)
            min_distance_arg = np.argmin(distances_from_root)

            if i == min_distance_arg:
                out.append(root)
                distances.append(distances_from_root[min_distance_arg])
            # The values in distances[] are a diagnostic on how good the
            # numerical accuracy is.

        # If the lens equation is solved correctly, there should be
        # either 3 or 5 solutions (corresponding to 3 or 5 images)
        if len(out) not in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            msg = ("Wrong number of solutions to the lens equation of binary" +
                   " lens.\nGot {:} and expected 3 or 5.\nThe parameters " +
                   "(m1, m2, s, source_x, source_y, solver) are:\n" +
                   "{:} {:} {:} {:} {:}  {:}\n\n" +
                   "Consider using 'point_source_point_lens' method for " +
                   "epochs when the source is very far from the lens. Note " +
                   "that it's different from 'point_source' method.")
            txt = msg.format(
                len(out), repr(self.mass_1), repr(self.mass_2),
                repr(self.separation), repr(source_x), repr(source_y),
                self._solver)

            if self._solver != "Skowron_and_Gould_12":
                txt += (
                    "\n\nYou should switch to using Skowron_and_Gould_12" +
                    " polynomial root solver. It is much more accurate than " +
                    "numpy.polynomial.polynomial.polyroots(). " +
                    "Skowron_and_Gould_12 method is selected in automated " +
                    "way if VBBL is imported properly.")
            distance = sqrt(source_x**2 + source_y**2)
            if (self.mass_2 > 1.e-6 * self.mass_1 and
                    (distance < 15. or distance < 2. * self.separation)):
                txt += ("\n\nThis is surprising error - please contact code " +
                        "authors and provide the above error message.")
            elif distance > 200.:
                txt += ("\n\nYou try to calculate magnification at huge " +
                        "distance from the source and this is causing an " +
                        "error.")
            txt += "\nMulensModel version: {:}".format(mm.__version__)

            raise ValueError(txt)
        if return_distances:
            return (np.array(out), np.array(distances))
        else:
            return np.array(out)

    def _get_jacobian_determinant(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._verify_polynomial_roots(
            source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        add_1 = self.mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self.mass_2 / (self._position_z2 - roots_ok_bar)**2
        derivative = add_1 + add_2 - self.shear_G

        return (1. - self.convergence_K)**2 - (derivative *
                                               np.conjugate(derivative))

    def point_source_magnification(self, source_x, source_y, vbbl_on=True):
        """
        Calculate point source magnification for given position. The
        origin of the coordinate system is at the center of mass and
        both masses are on X axis with higher mass at negative X; this
        means that the higher mass is at (X, Y)=(-s*q/(1+q), 0) and
        the lower mass is at (s/(1+q), 0).

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

        Returns :
            magnification: *float*
                Point source magnification.
        """
        if self._use_planet_frame:
            x_shift = -self.mass_1 / (self.mass_1 + self.mass_2)
        else:
            x_shift = self.mass_2 / (self.mass_1 + self.mass_2) - 0.5
        x_shift *= self.separation
        # We need to add this because in order to shift to correct frame.

        # Run point source using faster VBBL method.
        if _vbbl_wrapped and vbbl_on:
            s = float(self.separation)
            q = float(self.mass_2 / self.mass_1)
            x = float(source_x)
            y = float(source_y)

            magnification = _vbbl_binary_mag_0(
                s, q, x, y, self.convergence_K, self.shear_G.real,
                self.shear_G.imag)
        else:
            # Run point source using slower numpy method.
            magnification = self._get_point_source_Witt_Mao_95(
                source_x=float(source_x)+x_shift, source_y=float(source_y))

        if magnification < 1.:
            msg = (
                "error in BinaryLensWithShear.point_source_magnification()\n"
                "input:\n")
            params = [s, q, x, y, self.convergence_K, self.shear_G.real,
                      self.shear_G.imag, vbbl_on, _vbbl_wrapped]
            msg += " ".join([str(p) for p in params])
            msg += "\noutput: {:}".format(magnification)
            raise ValueError(msg)

        return magnification
