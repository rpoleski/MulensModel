import numpy as np
from astropy import __version__ as astropy__version__

MAG_ZEROPOINT = 22. # Defines magnitude at which flux = 1.

class Utils(object):
    '''a number of small functions used in different places'''

    def get_flux_from_mag(mag):
        '''transform magnitudes into fluxes'''
        flux = 10. ** (0.4 * (MAG_ZEROPOINT - mag))
        return flux
    get_flux_from_mag = staticmethod(get_flux_from_mag)

    def get_flux_and_err_from_mag(mag, err_mag):
        '''transform magnitudes into fluxes including errorbars'''
        flux = 10. ** (0.4 * (MAG_ZEROPOINT - mag))
        err_flux = err_mag * flux * np.log(10.) * 0.4
        return (flux, err_flux)
    get_flux_and_err_from_mag = staticmethod(get_flux_and_err_from_mag)

    def get_mag_from_flux(flux):
        '''transform fluxes into magnitudes'''
        mag = MAG_ZEROPOINT - 2.5 * np.log10(flux)
        return mag
    get_mag_from_flux = staticmethod(get_mag_from_flux)

    def get_mag_and_err_from_flux(flux, err_flux):
        '''transform fluxes into magnitudes including errorbars'''
        mag = MAG_ZEROPOINT - 2.5 * np.log10(flux)
        err_mag = (err_flux / flux) * 2.5 / np.log(10.)
        return (mag, err_mag)
    get_mag_and_err_from_flux = staticmethod(get_mag_and_err_from_flux)

    def astropy_version_check(minimum):
        '''check if astropy is installed at given or later version (input as a string)'''
        current = astropy__version__.split(".")
        required = minimum.split(".")
        for i in range(len(required)):
            if int(current[i]) < int(required[i]):
                return False
        return True

