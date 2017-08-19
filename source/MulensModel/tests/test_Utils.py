import numpy as np
import os

import MulensModel
import MulensModel.utils


SAMPLE_FILE_01 = os.path.join(MulensModel.MODULE_PATH, 
                                os.path.join("data", "phot_ob08092_O4.dat"))

def test_complex_fsum_1():
    z = [(0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1-1e+99j), (0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1+0.1j), (0.1+1e+99j), (0.1+0.1j)]
    assert MulensModel.utils.Utils.complex_fsum(z) == (1 + 0.8j)

def do_mag2flux_conversions_test(mag, mag_err):
    (flux, flux_err) = MulensModel.utils.Utils.get_flux_and_err_from_mag(
        mag, mag_err)
    (new_mag, new_mag_err) = MulensModel.utils.Utils.get_mag_and_err_from_flux(
        flux, flux_err)
    np.testing.assert_almost_equal(new_mag, mag, decimal=6)
    np.testing.assert_almost_equal(new_mag_err, mag_err, decimal=6)
    
def test_mag_and_flux_conversions_1():
    mag = 22.
    mag_err = 0.01

    (flux, flux_err) = MulensModel.utils.Utils.get_flux_and_err_from_mag(
        mag, mag_err)

    assert flux == 1.
    do_mag2flux_conversions_test(mag, mag_err)

def test_mag_and_flux_conversions_2():
    for mag in np.arange(22., 15., -1.):
        for mag_err in [0.01, 0.001, 0.1]:
            do_mag2flux_conversions_test(mag,mag_err)

