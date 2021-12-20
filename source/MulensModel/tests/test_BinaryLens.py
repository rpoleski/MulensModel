import sys
import numpy as np
import unittest

import MulensModel as mm


def test_small_q():
    """
    run calculations with small mass-ratio
    """
    q = 1.e-8
    s = 1.8
    x = s - 1. / s
    y = 5.e-7 * 0.

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    lens = mm.BinaryLens(m_1, m_2, s)
    result = lens.point_source_magnification(x, y)
    np.testing.assert_almost_equal(result, 3.6868957, decimal=3)


def test_vbbl_vs_point_source():
    """
    check if vbbl and point source calculations give the same result
    """
    s = 1.1
    q = 0.2
    x = s - 1. / s
    y = 5.e-7 * 0.

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    lens = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.1, shear_G=complex(0.1, -0.1))
    result = lens.point_source_magnification(x, y)
    result_vbbl = lens.vbbl_magnification(x, y, 1e-5)
    np.testing.assert_almost_equal(result, result_vbbl, decimal=3)


def test_standard_vs_shear():
    """
    check if standard and 0 shear calculations give the same result
    """
    s = 1.8
    q = 1e-6
    x = s - 1. / s
    y = 5.e-7 * 0.

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    lens = mm.BinaryLensWithShear(m_1, m_2, s,
                                  convergence_K=0.0, shear_G=complex(0, 0))
    lens_standard = mm.BinaryLens(m_1, m_2, s)
    result = lens.point_source_magnification(x, y)
    result_standard = lens_standard.point_source_magnification(x, y)
    np.testing.assert_almost_equal(result, result_standard, decimal=3)
# For q < 1e-7 the standard and shear methods start to deviate.


def test_binary_lens_hexadecapole():
    """
    tests hexadecapole and quadrupole calculation for planetary case
    with rho=0.001 and 3 values of gamma: 0, 0.5, and 1.0
    """
    s = 1.35
    q = 0.00578
    rho = 0.001
    x = 0.7142010570568691 - s * q / (1. + q)
    y = 0.00189679191923936

    pspl_mag = 4.691830779895085
    reference_00 = [5.017252440557196, 4.9587638949353146, pspl_mag]
    reference_05 = [4.981368071884021, 4.932070583431292, pspl_mag]
    reference_10 = [4.9454837032108445, 4.905377271927269, pspl_mag]

    m1 = 1. / (1. + q)
    m2 = q / (1. + q)

    bl = mm.BinaryLens(m1, m2, s)

    result_00 = bl.hexadecapole_magnification(
        x, y, rho=rho, gamma=0.0, all_approximations=True)
    result_05 = bl.hexadecapole_magnification(
        x, y, rho=rho, gamma=0.5, all_approximations=True)
    result_10 = bl.hexadecapole_magnification(
        x, y, rho=rho, gamma=1.0, all_approximations=True)

    np.testing.assert_almost_equal(result_00, reference_00)
    np.testing.assert_almost_equal(result_05, reference_05)
    np.testing.assert_almost_equal(result_10, reference_10)


def test_vbbl_1():
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLens(m_1, m_2, s)

    result = bl.vbbl_magnification(0.01, 0.01, 0.01)
    np.testing.assert_almost_equal(result, 18.2834436, decimal=3)


def test_vbbl_shear_1():
    """
    shear always uses a point source, so this test is not really fair
    """
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.0, shear_G=complex(0, 0))

    result = bl.vbbl_magnification(0.01, 0.01, 0.01)
    np.testing.assert_almost_equal(result, 18.2834436, decimal=1)


def test_ac_1():
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLens(m_1, m_2, s)

    result = bl.adaptive_contouring_magnification(
        0.01, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
    np.testing.assert_almost_equal(result, 18.2834436, decimal=3)
