import numpy as np

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
    """
    check basic magnification calculation using VBBL
    """
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLens(m_1, m_2, s)

    result = bl.vbbl_magnification(0.01, 0.01, 0.01)
    np.testing.assert_almost_equal(result, 18.2834436, decimal=3)


def test_vbbl_2():
    """
    Check VBBL magnification calculation for binary lens with finite source
    that was producing wrong results for VBBL earlier than v3.5
    """
    x = [-2.8798499936424813, -2.87980198609534, -2.879750341503788]
    y = [0.2603315602357186, 0.26034667859291694, 0.26036294250727565]
    s = 0.3121409537799967
    q = 0.0018654668855723224
    rho = 0.002966662955047919

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLens(m_1, m_2, s)
    results = [bl.vbbl_magnification(x_, y_, rho) for (x_, y_) in zip(x, y)]
    # VBBL 2.0.1 was returning:
    # [1.345365452870409, 1.368843518228974, 1.3442156685350604]
    # i.e., the second value was wrong.
    # VBBL 3.5 returns:
    # [1.3455151798453464, 1.3449680490180915, 1.344226470628085]
    np.testing.assert_almost_equal(results[2], 1.344226470628085, decimal=5)
    mean = (results[0] + results[2]) / 2
    np.testing.assert_almost_equal(results[1], mean, decimal=4)


def test_ac():
    """
    check basic magnification calculation using AdaptiveContouring
    """
    s = 0.8
    q = 0.1

    m_1 = 1. / (1. + q)
    m_2 = q / (1. + q)
    bl = mm.BinaryLens(m_1, m_2, s)

    # Previous version was failing because of bug in AC:
    # result = bl.adaptive_contouring_magnification(
    #     0.01, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
    # np.testing.assert_almost_equal(result, 18.2834436, decimal=3)
    result = bl.adaptive_contouring_magnification(
        0.06, 0.01, 0.01, accuracy=0.019, ld_accuracy=1e-3)
    np.testing.assert_almost_equal(result, 11.403036510555962, decimal=3)
