import numpy as np
import os

import MulensModel as mm

from test_Caustics import get_index_of_nearest_point

FILE_02 = os.path.join(mm.DATA_PATH, "unit_test_files",
                       "caustic_point_lens_with_shear.dat")
test_caustics_point = np.genfromtxt(FILE_02, names=['X', 'Y'], dtype=None)


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
