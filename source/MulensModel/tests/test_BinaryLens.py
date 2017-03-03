import sys
import numpy as np
import unittest

import MulensModel


def test_binary_lens_hexadecapole():
    """
    tests hexadecapole and quadrupole calculation for planetary case 
    with rho=0.001 and 3 values of gamma: 0, 0.5, and 1.0
    """
    s = 1.35
    q = 0.00578
    rho = 0.001
    x = 0.7142010570568691 
    y = 0.00189679191923936

    pspl_mag = 4.691830779895085
    reference_00 = [5.017252440557196, 4.9587638949353146, pspl_mag]
    reference_05 = [4.981368071884021, 4.932070583431292, pspl_mag]
    reference_10 = [4.9454837032108445, 4.905377271927269, pspl_mag]

    m1 = 1. / (1. + q)
    m2 = q / (1. + q)

    bl = MulensModel.BinaryLens(m1, m2, s)

    result_00 = bl.hexadecapole_magnification(x, y, rho=rho, 
                                        gamma=0.0, all_approximations=True)
    result_05 = bl.hexadecapole_magnification(x, y, rho=rho,
                                        gamma=0.5, all_approximations=True)
    result_10 = bl.hexadecapole_magnification(x, y, rho=rho,
                                        gamma=1.0, all_approximations=True)

    np.testing.assert_almost_equal(result_00, reference_00)
    np.testing.assert_almost_equal(result_05, reference_05)
    np.testing.assert_almost_equal(result_10, reference_10)

