import numpy as np

import MulensModel as mm


def test_ps_shear_1():
    """
    Test if vbbl_magnification with shear 0 and convergence 0 gives the same result 
    as point_source_magnification
    """
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.0, shear_G=complex(0, 0))

    result = bl.point_source_magnification(0.01, 0.01)
    np.testing.assert_almost_equal(result/18.2834436, 1., decimal=2)


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
