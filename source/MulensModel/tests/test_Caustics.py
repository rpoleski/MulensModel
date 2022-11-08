import numpy as np
import os

import MulensModel as mm


dir_ = os.path.join(mm.DATA_PATH, "unit_test_files")
FILE_01 = os.path.join(dir_, "MB11293_caustics.dat")
FILE_02 = os.path.join(dir_, "caustic_point_lens_with_shear.dat")

test_caustics = np.genfromtxt(FILE_01, names=['X', 'Y'], dtype=None)
test_caustics_point = np.genfromtxt(FILE_02, names=['X', 'Y'], dtype=None)


def get_index_of_nearest_point(test_data, x, y):
    """
    Find index of test_data (*dict*) that is closes to (x, y)
    """
    return np.argmin((test_data['X']-x)**2+(test_data['Y']-y)**2)


def test_caustic():
    """
    Make sure the caustic is properly calculated.
    """
    s = 0.548
    q = 0.0053

    caustics = mm.Caustics(q=q, s=s)

    x, y = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = get_index_of_nearest_point(test_caustics, x[i], y[i])

        np.testing.assert_almost_equal(
            x[i], test_caustics['X'][index], decimal=5)
        np.testing.assert_almost_equal(
            y[i], test_caustics['Y'][index], decimal=5)


def test_point_lens_caustic():
    """
    Test point lens caustics with external mass sheet
    """
    convergence_K = 0.1
    shear_G = complex(-0.1, -0.2)

    caustics = mm.CausticsPointWithShear(convergence_K=convergence_K,
                                         shear_G=shear_G)

    (x, y) = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = get_index_of_nearest_point(test_caustics_point, x[i], y[i])

        np.testing.assert_almost_equal(x[i], test_caustics_point['X'][index])
        np.testing.assert_almost_equal(y[i], test_caustics_point['Y'][index])


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
