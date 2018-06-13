import sys
import os
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

import MulensModel
from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit


def test_fit_get_input_format():
    """read sample file and get brightness in its original format"""
    n = 100
    mag = 15.
    time = np.ones(n) * 2456789.
    dataset = MulensData(data_list=[time, time*0.+mag, time*0.+0.01])
    fit = Fit(data=[dataset], magnification=[np.ones(n)])
    fit.fit_fluxes()

    input_fmt = fit.get_input_format(data=dataset)
    err_msg = 'Fit.get_input_format() returns something wrong'
    np.testing.assert_almost_equal(input_fmt, mag, err_msg=err_msg)
   
