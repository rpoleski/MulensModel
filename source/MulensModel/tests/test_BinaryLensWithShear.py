import numpy as np

import MulensModel as mm


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
    np.testing.assert_almost_equal(result/18.2834436, 1., decimal=2)


# Test binary lens with shear magnification against ray tracing in
# Vedantham et al. 2017 (ApJ 845 2 89).
def test_vbbl_shear_2():
    """
    shear always uses a point source
    """
    s = 0.8
    q = 0.9

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.05, shear_G=complex(0.1, -0.06))

    y = 0.0
    x = 0.11776447105788423
    result = bl.vbbl_magnification(x, y, 0.0001)
    np.testing.assert_almost_equal(result, 15.07985008909832)


def test_vbbl_shear_3():
    """
    shear always uses a point source
    """
    s = 0.7
    q = 0.01

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLensWithShear(
        m_1, m_2, s, convergence_K=0.03, shear_G=complex(-0.07, 0.03))

    y = 0.2
    x = 0.03792415169660668
    result = bl.vbbl_magnification(x, y, 0.0001)
    np.testing.assert_almost_equal(result, 6.480196704201193)
