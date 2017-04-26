import numpy as np

import MulensModel
from MulensModel.caustics import Caustics


MODULE_PATH = "/".join(MulensModel.__file__.split("/source")[:-1])

SAMPLE_FILE_01 = MODULE_PATH + "/data/MB11293_caustics.dat"

test_caustics = np.genfromtxt(SAMPLE_FILE_01, names=['X', 'Y'], dtype=None)


def test_caustic():
    s = 0.548
    q = 0.0053

    caustics = Caustics(q=q, s=s)
    caustics.calculate()

    shift_x = s / 2.
    
    np.testing.assert_almost_equal(caustics._x[100], test_caustics['X'][100]+shift_x, decimal=3)
    np.testing.assert_almost_equal(caustics._y[100], test_caustics['Y'][100], decimal=3)

