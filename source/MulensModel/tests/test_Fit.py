import sys
import os
import numpy as np
from numpy.testing import assert_almost_equal as almost

from astropy.coordinates import SkyCoord
from astropy import units as u

import MulensModel as mm


def test_fit_get_input_format():
    """read sample file and get brightness in its original format"""
    n = 100
    mag = 15.
    time = np.ones(n) * 2456789.
    dataset = mm.MulensData(data_list=[time, time*0.+mag, time*0.+0.01],
                            phot_fmt='flux')
    fit = mm.Fit(data=[dataset], magnification=[np.ones(n)])
    fit.fit_fluxes()

    input_fmt = fit.get_input_format(data=dataset)
    err_msg = 'Fit.get_input_format() returns something wrong'
    almost(input_fmt, mag, err_msg=err_msg)


def test_update():
    """
    check that Fit.update() works properly
    """
    n = 100
    (flux_0, flux_1) = (10., 20.)
    time = np.linspace(2454000., 2454100., n)
    magnification = (time - 2454000.) / 100.
    flux = flux_0 + flux_1 * magnification
    err = flux * 0. + 1.e-3
    dataset_1 = mm.MulensData(data_list=[time, flux, err], phot_fmt='flux')
    fit_1 = mm.Fit(data=[dataset_1], magnification=magnification)
    fit_1.fit_fluxes()

    n = 200
    (flux_2, flux_3) = (15., 2.)
    magnification = (time - 2454050) / 30.
    flux = flux_2 + flux_3 * magnification
    err = flux * 0. + 1.e-3
    dataset_2 = mm.MulensData(data_list=[time, flux, err], phot_fmt='flux')
    fit_2 = mm.Fit(data=[dataset_2], magnification=magnification)
    fit_2.fit_fluxes()

    fit_1.update(fit_2)

    almost(fit_1.flux_of_sources(dataset_1), flux_1)
    almost(fit_1.flux_of_sources(dataset_2), flux_3)
    almost(fit_2.flux_of_sources(dataset_2), flux_3)

    almost(fit_1.blending_flux(dataset_1), flux_0)
    almost(fit_1.blending_flux(dataset_2), flux_2)
    almost(fit_2.blending_flux(dataset_2), flux_2)
