import numpy as np

import MulensModel as mm


def test_ps_shear_1():
    """
    Test if vbbl_magnification() with shear 0 and convergence 0
    gives the same result as point_source_magnification()
    """
    s = 0.8
    q = 0.1
    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    blws = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.0, shear_G=complex(0, 0))
    result = blws.point_source_magnification(0.01, 0.01)

    bl = mm.BinaryLens(m_1, m_2, s)
    result_ref = bl.point_source_magnification(0.01, 0.01)

    np.testing.assert_almost_equal(result, result_ref)


def test_ps_shear_2():
    """
    Test binary lens with shear magnification against ray tracing in
    Vedantham et al. 2017 (ApJ 845 2 89). Two different tests of
    parameters are used below.
    """
    s = 0.8
    q = 0.9
    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.05, shear_G=complex(0.1, -0.06))

    result = bl.point_source_magnification(0.11776447105788423, 0.)
    np.testing.assert_almost_equal(result, 15.07985008909832)

    s = 0.7
    q = 0.01
    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.03, shear_G=complex(-0.07, 0.03))

    result = bl.point_source_magnification(0.03792415169660668, 0.2)
    np.testing.assert_almost_equal(result, 6.480196704201193)


def test_vbbl_vs_MM():
    """
    check if implementations in VBBL and MM give the same result
    """
    s = 1.1
    q = 0.2
    x = s - 1. / s
    y = 0.

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    lens = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.1, shear_G=complex(0.1, -0.1))
    result = lens.point_source_magnification(x, y, vbbl_on=False)
    result_vbbl = lens.point_source_magnification(x, y, vbbl_on=True)
    np.testing.assert_almost_equal(result, result_vbbl)


def test_standard_vs_shear():
    """
    check if standard and 0 shear calculations give the same result
    """
    s = 1.8
    q = 1e-6
    x = s - 1. / s
    y = 0.

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)

    lens = mm.BinaryLensWithShear(m_1, m_2, s,
                                  convergence_K=0.0, shear_G=complex(0, 0))
    lens_standard = mm.BinaryLens(m_1, m_2, s)
    result = lens.point_source_magnification(x, y)
    result_standard = lens_standard.point_source_magnification(x, y)
    np.testing.assert_almost_equal(result, result_standard, decimal=5)
