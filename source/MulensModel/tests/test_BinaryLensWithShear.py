import numpy as np

import MulensModel as mm
from test_BinaryLens import make_trajectory


def test_ps_shear_1():
    """
    Test if vbbl_magnification() with shear 0 and convergence 0
    gives the same result as point_source_magnification()
    """
    s = 0.8
    q = 0.1
    x, y = 0.01, 0.01

    trajectory = make_trajectory((x, y), {'s': s, 'q': q})

    blws = mm.BinaryLensPointSourceWithShearVBBLMagnification(
        trajectory=trajectory, convergence_K=0.0, shear_G=complex(0, 0))
    result = blws.get_magnification()

    bl = mm.BinaryLensPointSourceWM95Magnification(
        trajectory=trajectory)
    result_ref = bl.get_magnification()

    np.testing.assert_almost_equal(result, result_ref)


def test_ps_shear_2():
    """
    Test binary lens with shear magnification against ray tracing in
    Vedantham et al. 2017 (ApJ 845 2 89). Two different tests of
    parameters are used below.
    """
    s = 0.8
    q = 0.9
    x, y = 0.11776447105788423, 0.

    trajectory = make_trajectory((x, y), {'s': s, 'q': q})

    expected = 15.07985008909832
    bl = mm.BinaryLensPointSourceWithShearVBBLMagnification(
        trajectory=trajectory, convergence_K=0.05,
        shear_G=complex(0.1, -0.06))
    bl_wm95 = mm.BinaryLensPointSourceWithShearWM95Magnification(
        trajectory=trajectory, convergence_K=0.05,
        shear_G=complex(0.1, -0.06))

    result = bl.get_magnification()
    np.testing.assert_almost_equal(result, expected)
    np.testing.assert_almost_equal(bl_wm95.get_magnification(), expected)


def test_ps_shear_3():
    s = 0.7
    q = 0.01
    x, y = 0.03792415169660668, 0.2
    expected = 6.480196704201193

    trajectory = make_trajectory((x, y), {'s': s, 'q': q})

    bl = mm.BinaryLensPointSourceWithShearVBBLMagnification(
        trajectory=trajectory, convergence_K=0.03, shear_G=complex(-0.07, 0.03))
    bl_wm95 = mm.BinaryLensPointSourceWithShearWM95Magnification(
        trajectory=trajectory, convergence_K=0.03, shear_G=complex(-0.07, 0.03))

    result = bl.get_magnification()
    np.testing.assert_almost_equal(result, expected)
    np.testing.assert_almost_equal(bl_wm95.get_magnification(), expected)


def test_vbbl_vs_MM():
    """
    check if implementations in VBBL and MM give the same result
    """
    s = 1.1
    q = 0.2
    x = s - 1. / s
    y = 0.

    trajectory = make_trajectory((x, y), {'s': s, 'q': q})

    lens_wm95 = mm.BinaryLensPointSourceWithShearVBBLMagnification(
       trajectory=trajectory, convergence_K=0.1, shear_G=complex(0.1, -0.1))
    result = lens_wm95.get_magnification()

    lens = mm.BinaryLensPointSourceWithShearWM95Magnification(
       trajectory=trajectory, convergence_K=0.1, shear_G=complex(0.1, -0.1))
    result_vbbl = lens.get_magnification()

    np.testing.assert_almost_equal(result, result_vbbl)


def test_standard_vs_shear():
    """
    check if standard and 0 shear calculations give the same result
    """
    s = 1.8
    q = 1e-6
    x = s - 1. / s
    y = 0.

    trajectory = make_trajectory((x, y), {'s': s, 'q': q})

    lens = mm.BinaryLensPointSourceWithShearVBBLMagnification(
        trajectory=trajectory, convergence_K=0.0, shear_G=complex(0, 0))
    lens_standard = mm.BinaryLensPointSourceWM95Magnification(
        trajectory=trajectory)

    result = lens.get_magnification()
    result_standard = lens_standard.get_magnification()
    np.testing.assert_almost_equal(result, result_standard, decimal=5)
