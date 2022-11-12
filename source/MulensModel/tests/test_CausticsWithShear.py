import numpy as np
import os

import MulensModel as mm


SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "unit_test_files", "MB11293_caustics.dat")

test_caustics = np.genfromtxt(SAMPLE_FILE_01, names=['X', 'Y'], dtype=None)


def test_caustic_shear():
    """
    Trivial test of shear and convergence equal to 0 - same as for Caustics.
    """
    s = 0.548
    q = 0.0053

    caustics = mm.CausticsWithShear(
        q=q, s=s, convergence_K=0.0, shear_G=complex(0, 0))

    x, y = caustics.get_caustics(n_points=1000)
    for i in range(0, len(x), 100):
        index = np.argmin(
            np.sqrt((test_caustics['X']-x[i])**2+(test_caustics['Y']-y[i])**2))

        np.testing.assert_almost_equal(
            x[i], test_caustics['X'][index], decimal=5)
        np.testing.assert_almost_equal(
            y[i], test_caustics['Y'][index], decimal=5)
