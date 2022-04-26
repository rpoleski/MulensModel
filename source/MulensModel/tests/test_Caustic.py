import numpy as np
import os

import MulensModel as mm
from MulensModel.caustics import Caustics


SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "unit_test_files", "MB11293_caustics.dat")

test_caustics = np.genfromtxt(SAMPLE_FILE_01, names=['X', 'Y'], dtype=None)


def test_caustic():
    s = 0.548
    q = 0.0053

    caustics = mm.Caustics(q=q, s=s)

    x, y = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = np.argmin(
            np.sqrt((test_caustics['X']-x[i])**2+(test_caustics['Y']-y[i])**2))

        np.testing.assert_almost_equal(
            x[i], test_caustics['X'][index], decimal=5)
        np.testing.assert_almost_equal(
            y[i], test_caustics['Y'][index], decimal=5)

# def test_point_caustic():
#     """
#     Test point lens caustics with external mass sheet
#     """
#     convergence_K = 0.1
#     shear_G = complex(-0.1,-0.2)

#     caustics = mm.CausticsPointWithShear(convergence_K=convergence_K, shear_G=shear_G)

#     x, y = caustics.get_caustics(n_points=1000)
#     for i in range(0, len(x), 100):
#         index = np.argmin(
#             np.sqrt((test_caustics['X']-x[i])**2+(test_caustics['Y']-y[i])**2))

#         np.testing.assert_almost_equal(
#             x[i], test_caustics['X'][index], decimal=5)
#         np.testing.assert_almost_equal(
#             y[i], test_caustics['Y'][index], decimal=5)

def test_point_caustic_no_external_mass_sheet():
    """
    Test point lens caustics with external mass sheet reduces
    to no external mass sheet when convergence and shear are zero
    """
    convergence_K = 0.0
    shear_G = complex(0.0, 0.0)

    caustics = mm.CausticsPointWithShear(convergence_K=convergence_K,
                                         shear_G=shear_G)
    (x, y) = caustics.get_caustics(n_points=100)
    for i in range(0, len(x), 100):
        np.testing.assert_almost_equal(x[i], 0.0, decimal=5)
        np.testing.assert_almost_equal(y[i], 0.0, decimal=5)
