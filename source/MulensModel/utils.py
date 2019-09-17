import numpy as np
from math import fsum, pow, sqrt
import warnings

from astropy import __version__ as astropy__version__
from astropy.time import Time
from astropy.coordinates.builtin_frames.utils import get_jd12
from astropy import _erfa as erfa


MAG_ZEROPOINT = 22.  # Defines magnitude at which flux = 1.

month_3letter_to_2digit = {
    'Jan': '01',
    'Feb': '02',
    'Mar': '03',
    'Apr': '04',
    'May': '05',
    'Jun': '06',
    'Jul': '07',
    'Aug': '08',
    'Sep': '09',
    'Oct': '10',
    'Nov': '11',
    'Dec': '12'
    }


class Utils(object):
    """ A number of small functions used in different places """

    def get_flux_from_mag(mag, zeropoint=MAG_ZEROPOINT):
        """ transform magnitudes into fluxes """
        flux = 10. ** (0.4 * (zeropoint - mag))
        return flux
    get_flux_from_mag = staticmethod(get_flux_from_mag)

    def get_flux_and_err_from_mag(mag, err_mag, zeropoint=MAG_ZEROPOINT):
        """ transform magnitudes into fluxes including errorbars """
        flux = 10. ** (0.4 * (zeropoint - mag))
        err_flux = err_mag * flux * np.log(10.) * 0.4
        return (flux, err_flux)
    get_flux_and_err_from_mag = staticmethod(get_flux_and_err_from_mag)

    def get_mag_from_flux(flux, zeropoint=MAG_ZEROPOINT):
        """ transform fluxes into magnitudes """
        if np.any(flux <= 0.):
            warnings.warn(
                "Flux to magnitude conversion approached negative flux",
                UserWarning)
        mag = zeropoint - 2.5 * np.log10(flux)
        return mag
    get_mag_from_flux = staticmethod(get_mag_from_flux)

    def get_mag_and_err_from_flux(flux, err_flux, zeropoint=MAG_ZEROPOINT):
        """ transform fluxes into magnitudes including errorbars """
        if np.any(flux <= 0.):
            warnings.warn(
                "Flux to magnitude conversion approached negative flux",
                UserWarning)
        mag = zeropoint - 2.5 * np.log10(flux)
        err_mag = (err_flux / flux) * 2.5 / np.log(10.)
        return (mag, err_mag)
    get_mag_and_err_from_flux = staticmethod(get_mag_and_err_from_flux)

    def complex_fsum(arguments):
        """
        accurate floating points sum of complex numbers in iterable
        arguments
        """
        real = [arg.real for arg in arguments]
        imag = [arg.imag for arg in arguments]
        return fsum(real) + fsum(imag) * 1j
    complex_fsum = staticmethod(complex_fsum)

    def dot(cartesian, vector):
        """ dot product of Astropy CartersianRepresentation and np.ndarray """
        return (cartesian.x * vector[0] + cartesian.y * vector[1] +
                cartesian.z * vector[2])
    dot = staticmethod(dot)

    def vector_product_normalized(vector_1, vector_2):
        """
        get vector that is perpendicular to the 2 above and is
        normalized
        """
        vector_product = np.cross(vector_1, vector_2)
        return vector_product / np.linalg.norm(vector_product)
    vector_product_normalized = staticmethod(vector_product_normalized)

    # The two functions below implement convention introduced by:
    # An et al. 2002 (ApJ 572, 521)
    # https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract
    def gamma_to_u(gamma):
        """ transform gamma limb darkening coefficient to u """
        return 3. * gamma / (2. + gamma)
    gamma_to_u = staticmethod(gamma_to_u)

    def u_to_gamma(u):
        """ transform u limb darkening coefficient to gamma """
        return (2. * u) / (3. - u)
    u_to_gamma = staticmethod(u_to_gamma)

    def get_n_caustics(s, q):
        """
        Find number of caustics for binary lens.

        Parameters :
            s: *float*
                Separation of binary lens relative to theta_E.

            q: *float*
                Mass ratio (<1).

        Returns :
            n_caustics: *int*
                Number of caustics: 1, 2, or 3.
        """
        limit = (1. + q) / (1. + q**(1./3.))**3
        if s > 1. / sqrt(limit):
            return 2
        elif s < pow(limit, 0.25):
            return 3
        else:
            return 1

    def date_change(text):
        """
        changes format: '2015-Oct-30 12:00' -> '2015-10-30 12:00'
        """
        text = text.decode('UTF-8')
        str_components = text.split('-')
        if len(str_components) == 1:
            raise ValueError("Can't run date_change() for {:}".format(text))
        return '-'.join((
            str_components[0], month_3letter_to_2digit[str_components[1]],
            str_components[2]))
    date_change = staticmethod(date_change)

    def velocity_of_Earth(full_BJD):
        """
        Calculate 3D velocity of Earth for given epoch.

        If you need velocity projected on the plane of the sky, then use
        :py:func:`~MulensModel.coordinates.Coordinates.v_Earth_projected`

        Parameters :
            full_BJD: *float*
                Are we fitting for blending flux? If not then blending flux is
                fixed to 0.  Default is the same as
                :py:func:`MulensModel.fit.Fit.fit_fluxes()`.

        Returns :
            velocity: *np.ndarray* (*float*, size of (3,))
                3D velocity in km/s. The frame follows *Astropy* conventions.
        """
        # The 4 lines below, that calculate velocity for given epoch,
        # are based on astropy 1.3 code:
        # https://github.com/astropy/astropy/blob/master/astropy/
        # coordinates/solar_system.py
        time = Time(full_BJD, format='jd', scale='tdb')
        (jd1, jd2) = get_jd12(time, 'tdb')
        (earth_pv_helio, earth_pv_bary) = erfa.epv00(jd1, jd2)
        factor = 1731.45683  # This scales AU/day to km/s.
        # The returned values are of np.ndarray type in astropy v1 and v2,
        # but np.void in v3. The np.asarray() works in both cases.
        velocity = np.asarray(earth_pv_bary[1]) * factor
        return velocity
    velocity_of_Earth = staticmethod(velocity_of_Earth)

    def astropy_version_check(minimum):
        """
        check if astropy is installed at given or later version (input
        as a string)
        """
        current = astropy__version__.split(".")
        required = minimum.split(".")
        for i in range(len(required)):
            if int(current[i]) < int(required[i]):
                return False
        return True
