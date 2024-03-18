from numpy.testing import assert_almost_equal
import warnings

import MulensModel as mm


# The test presented check if VBBL code was imported properly and
# they intentionally call private functions from mm.binarylensimports.

# The values for tests were copied from other tests.


def test_VBBL_dark():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_binary_mag_dark()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_binary_mag_dark(
        1.35, 0.00578, 0.6598560303179819, -0.05280389642190758, 0.01,
        0.001, 0.0)
    assert_almost_equal(out, 1.635980468652538)


def test_VBBL_finite():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_binary_mag_finite()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_binary_mag_finite(
        0.8, 0.1, 0.01, 0.01, 0.01, 0.001)
    assert_almost_equal(out, 18.283392940574107)


def test_VBBL_point():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_binary_mag_point()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_binary_mag_point(0.8, 0.1, 0.01, 0.01)
    assert_almost_equal(out, 18.185448359975954)


def test_VBBL_shear():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_binary_mag_point_shear()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_binary_mag_point_shear(
        1.1, 0.2, 0.19090909090909103, 0.0, 0.1, 0.1, -0.1)
    assert_almost_equal(out, 7.797177952275903)


def test_VBBL_SG12_5():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_SG12_5()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_SG12_5(
        -0.17073496203535826, -0.2960054288981844, 0.2635469112719409,
        0.8257865631242759, 0.8806819361918513, 0.31179330108200043,
        0.029491785441276515, 0.17299410121242445, 0.2364363674161051,
        0.13879110194896355, 0.056792767235166075, -0.003964567527886535)
    expected = [
        -0.58000512, -1.6727006, -0.57653394, 0.56004431, -0.5526021,
        0.90541606, -0.16774112, -0.76401705, -0.14860692, -0.04307995]
    assert_almost_equal(out, expected)


def test_VBBL_SG12_9():
    """
    Directly (hence, calling private function) test imported VBBL function:
    _vbbl_SG12_9()
    """
    if not mm.binarylensimports._vbbl_wrapped:
        warnings.warn("VBBL not imported", UserWarning)
        return

    out = mm.binarylensimports._vbbl_SG12_9(
        -0.0006162037037037004, 0.018455861111111117, 0.06027712777777768,
        -0.08050824074074037, -0.30724011207070717, -0.5752337759963271,
        -0.524794951286975, -0.2288742865013774, -0.0350081818181818, 0.0,
        0.0006162037037037004, 0.005693722222222227, 0.006298813888888853,
        -0.0035301525925927266, 0.063350366010101, 0.1906118155555555,
        0.09094615694687937, 0.029791669421487667, 0.03019545454545458,
        0.0158)
    expected = [
        -1.50992662, 2.47560471, -1.61711916, -0.4453668, -0.75676062,
        -0.14298556, 0.36097464, -0.29933664, 0.02381135, -0.74906129,
        -3.63807783,  0.90118702, 1.29719497, 0.60196247, -0.63649524,
        0.0565556, -0.0138864,  -0.03508702]
    assert_almost_equal(out, expected)
