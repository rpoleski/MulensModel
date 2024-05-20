import numpy as np
import os

import MulensModel as mm


FILE_01 = os.path.join(mm.DATA_PATH, "unit_test_files",
                       "MB11293_caustics.dat")
test_caustics = np.genfromtxt(FILE_01, names=['X', 'Y'], dtype=None)


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

    caustics = mm.CausticsBinary(q=q, s=s)

    x, y = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = get_index_of_nearest_point(test_caustics, x[i], y[i])

        np.testing.assert_almost_equal(
            x[i], test_caustics['X'][index], decimal=5)
        np.testing.assert_almost_equal(
            y[i], test_caustics['Y'][index], decimal=5)
