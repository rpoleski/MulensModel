"""
File with general code used in other parts of MulensModel package.

Most importantly there are Utils and PlotUtils classes.
"""
import numpy as np
from math import fsum, pow, sqrt
import warnings
from matplotlib.colors import ColorConverter

from astropy import __version__ as astropy__version__
from astropy.time import Time
from astropy.coordinates.builtin_frames.utils import get_jd12
try:
    import erfa
except Exception:
    from astropy import _erfa as erfa
# In astropy 4.2 they removed astropy._erfa module and made a separate
# package erfa that is required by astropy.


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

    def get_flux_from_mag(mag, zeropoint=None):
        """
        Transform magnitudes into fluxes.

        Parameters :
            mag: *np.ndarray* or *float*
                Values to be transformed.

            zeropoint: *float*
                Zeropoint of magnitude scale.
                Defaults to 22. - double check if you want to change this.

        Returns :
            flux: *np.ndarray* or *float*
                Calculated fluxes. Type is the same as *mag* parameter.
        """
        if zeropoint is None:
            zeropoint = MAG_ZEROPOINT
        flux = 10. ** (0.4 * (zeropoint - mag))
        return flux
    get_flux_from_mag = staticmethod(get_flux_from_mag)

    def get_flux_and_err_from_mag(mag, err_mag, zeropoint=None):
        """
        Transform magnitudes and their uncertainties into flux space.

        Parameters :
            mag: *np.ndarray* or *float*
                Magnitude values to be transformed.

            err_mag: *np.ndarray* or *float*
                Uncertainties of magnitudes to be transformed.

            zeropoint: *float*
                Zeropoint of magnitude scale.
                Defaults to 22. - double check if you want to change this.

        Returns :
            flux: *np.ndarray* or *float*
                Calculated fluxes. Type is the same as *mag* parameter.

            err_flux: *np.ndarray* or *float*
                Calculated flux uncertainties. Type is *float* if both *mag*
                and *err_mag* are *floats* and *np.ndarray* otherwise.
        """
        if zeropoint is None:
            zeropoint = MAG_ZEROPOINT
        flux = 10. ** (0.4 * (zeropoint - mag))
        err_flux = err_mag * flux * np.log(10.) * 0.4
        return (flux, err_flux)
    get_flux_and_err_from_mag = staticmethod(get_flux_and_err_from_mag)

    def get_mag_from_flux(flux, zeropoint=None):
        """
        Transform fluxes into magnitudes.

        Parameters :
            flux: *np.ndarray* or *float*
                Values to be transformed.

            zeropoint: *float*
                Zeropoint of magnitude scale.
                Defaults to 22. - double check if you want to change this.

        Returns :
            mag: *np.ndarray* or *float*
                Calculated fluxes. Type is the same as *flux* parameter.
        """
        if zeropoint is None:
            zeropoint = MAG_ZEROPOINT
        if np.any(flux <= 0.):
            warnings.warn(
                "Flux to magnitude conversion approached negative flux",
                UserWarning)
        mag = zeropoint - 2.5 * np.log10(flux)
        return mag
    get_mag_from_flux = staticmethod(get_mag_from_flux)

    def get_mag_and_err_from_flux(flux, err_flux, zeropoint=None):
        """
        Transform fluxes and their uncertainties into magnitude space.

        Parameters :
            flux: *np.ndarray* or *float*
                Flux values to be transformed.

            err_flux: *np.ndarray* or *float*
                Uncertainties of fluxes to be transformed.

            zeropoint: *float*
                Zeropoint of magnitude scale.
                Defaults to 22. - double check if you want to change this.

        Returns :
            mag: *np.ndarray* or *float*
                Calculated fluxes. Type is the same as *mag* parameter.

            err_mag: *np.ndarray* or *float*
                Calculated magnitude uncertainties. Type is *float* if both
                *flux* and *err_flux* are *floats* and *np.ndarray* otherwise.
        """
        if zeropoint is None:
            zeropoint = MAG_ZEROPOINT
        if np.any(flux <= 0.):
            warnings.warn(
                "Flux to magnitude conversion approached negative flux",
                UserWarning)
        mag = zeropoint - 2.5 * np.log10(flux)
        err_mag = (err_flux / flux) * 2.5 / np.log(10.)
        return (mag, err_mag)
    get_mag_and_err_from_flux = staticmethod(get_mag_and_err_from_flux)

    # The two functions below implement convention introduced by:
    # An et al. 2002 (ApJ 572, 521)
    # https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract
    def gamma_to_u(gamma):
        """
        Transform gamma limb darkening coefficient to u in convention
        introduced by `An et al. 2008 (ApJ 681, 1593)
        <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_.

        Parameters :
            gamma: *float*
                Limb darkening coefficient in gamma convention.

        Returns :
            u: *float*
                Limb darkening coefficient in u convention.
        """
        return 3. * gamma / (2. + gamma)
    gamma_to_u = staticmethod(gamma_to_u)

    def u_to_gamma(u):
        """
        Transform u limb darkening coefficient to gamma in convention
        introduced by `An et al. 2008 (ApJ 681, 1593)
        <https://ui.adsabs.harvard.edu/abs/2002ApJ...572..521A/abstract>`_.

        Parameters :
            u: *float*
                Limb darkening coefficient in u convention.

        Returns :
            gamma: *float*
                Limb darkening coefficient in gamma convention.
        """
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
    get_n_caustics = staticmethod(get_n_caustics)

    def velocity_of_Earth(full_BJD):
        """
        Calculate 3D velocity of Earth for given epoch.

        If you need velocity projected on the plane of the sky, then use
        :py:func:`~MulensModel.coordinates.Coordinates.v_Earth_projected`

        Parameters :
            full_BJD: *float*
                Barycentric Julian Data. Full means it should begin
                with 245... or 246...

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

    def complex_fsum(arguments):
        """
        Accurate floating points sum of complex numbers in iterable.

        Parameters :
            arguments: *iterable* (e.g., *list* or *np.ndarray*)

        Returns :
            sum: *complex*
                Sum of *arguments*.
        """
        real = [arg.real for arg in arguments]
        imag = [arg.imag for arg in arguments]
        return fsum(real) + fsum(imag) * 1j
    complex_fsum = staticmethod(complex_fsum)

    def vector_product_normalized(vector_1, vector_2):
        """
        Get vector that is perpendicular to the 2 input ones and is normalized.

        Parameters :
            vector_1: *np.ndarray* (3,)
                First vector.

            vector_2: *np.ndarray* (3,)
                Second vector.

        Returns :
            vector_product_norm: *np.ndarray* (3,)
                Normalized vector product of the inputs.
        """
        vector_product = np.cross(vector_1, vector_2)
        return vector_product / np.linalg.norm(vector_product)
    vector_product_normalized = staticmethod(vector_product_normalized)

    def dot(cartesian, vector):
        """
        Dot product of 2 vectors represented by
        *astropy.CartesianRepresentation* and *np.ndarray*.

        Parameters :
            cartesian: *astropy.CartesianRepresentation*
                First vector.

            vector: *np.ndarray*
                Second vector (length 3).

        Returns :
            dot_product: *astropy.Quantity*
                Dot product of inputs.
        """
        return (cartesian.x * vector[0] + cartesian.y * vector[1] +
                cartesian.z * vector[2])
    dot = staticmethod(dot)

    def date_change(text):
        """
        Change format of month in date, e.g.
        *'2015-Oct-30 12:00'* -> *'2015-10-30 12:00'*

        Parameters :
            text: *str*
                Date to be changed.

        Returns :
            out_text: *str*
                Changed text.
        """
        text = text.decode('UTF-8')
        str_components = text.split('-')
        if len(str_components) == 1:
            raise ValueError("Can't run date_change() for {:}".format(text))
        return '-'.join((
            str_components[0], month_3letter_to_2digit[str_components[1]],
            str_components[2]))
    date_change = staticmethod(date_change)

    def astropy_version_check(minimum):
        """
        Check if *astropy* package is installed at given or later version.

        Parameters :
            minimum: *str*
                Minimum version, e.g., *'3.1.2'*.

        Returns :
            out: *bool*
                Is the installed version later?
        """
        current = astropy__version__.split(".")
        required = minimum.split(".")
        for i in range(len(required)):
            if int(current[i]) < int(required[i]):
                return False
        return True
    astropy_version_check = staticmethod(astropy_version_check)


class PlotUtils(object):
    """
    A number of small functions related to plotting used in different places
    """

    def get_y_value_y_err(phot_fmt, flux, flux_err):
        """
        Change the format of data for Y axis.

        Parameters :
            phot_fmt: *str*
                Requested format of output: 'mag' or 'flux'

            flux: *np.ndarray*
                fluxes values for Y axis

            flux_err: *np.ndarray*
                flux uncertainties for Y axis

        Returns :
            values: *np.ndarray*
                Values in the requested format

            uncertainties: *np.ndarray*
                Uncertainties in the requested format
        """
        if phot_fmt == 'mag':
            return Utils.get_mag_and_err_from_flux(flux, flux_err)
        elif phot_fmt == 'flux':
            return (flux, flux_err)
        else:
            raise ValueError(
                'Unrecognized photometry format: {:}, '.format(phot_fmt) +
                'allowed values are "mag" and "flux"')
    get_y_value_y_err = staticmethod(get_y_value_y_err)

    def find_subtract(subtract_2450000=False, subtract_2460000=False):
        """
        Find value that is supposed to be subtracted from time vector.

        Parameters :
            subtract_2450000: *bool*
                Should we subtract 2450000?

            subtract_2460000: *bool*
                Should we subtract 2460000?

        Returns :
            subtract: *float*
                Value to be subtracted.
        """
        subtract = 0.
        if subtract_2450000:
            if subtract_2460000:
                raise ValueError("subtract_2450000 and subtract_2460000 " +
                                 "cannot be both True")
            subtract = 2450000.
        if subtract_2460000:
            subtract = 2460000.

        return subtract
    find_subtract = staticmethod(find_subtract)

    def find_subtract_xlabel(subtract_2450000=False, subtract_2460000=False):
        """
        Find string for xlabel

        Parameters :
            subtract_2450000: *bool*
                Should we subtract 2450000?

            subtract_2460000: *bool*
                Should we subtract 2460000?

        Returns :
            xlabel: *str*
                String to be used by *pyplot.xlabel()*
        """
        if subtract_2450000:
            if subtract_2460000:
                raise ValueError('subtract_2450000 and subtract_2460000 ' +
                                 'cannot be both True')
            out = 'Time - 2450000'
        elif subtract_2460000:
            out = 'Time - 2460000'
        else:
            out = 'Time'

        return out
    find_subtract_xlabel = staticmethod(find_subtract_xlabel)

    def get_color_differences(color_list, color):
        # I couldn't make the link below to work in sphinx because of "_" - RP
        """
        Calculate color difference between a list of colors and a single color.
        Uses algorithm from Wikipedia page:
        https://en.wikipedia.org/wiki/Color_difference

        Parameters :
            color_list: *list* of *str*
                list of matplotlib colors e.g., ``['black', '#E9967A']``

            color: *str*
                single matplotlib color

        Returns :
            differences: *np.ndarray*
                differences of colors, values < 0.3 are very similar
        """
        rgba = ColorConverter.to_rgba
        array = np.array(
            [[float(x) for x in list(rgba(c))[:3]] for c in color_list])
        # We use float above because some versions of matplotlib return str.
        color_value = [float(x) for x in list(rgba(color))[:3]]
        mean_red = 0.5 * (array[:, 0] + color_value[0])
        diffs = (array - color_value)**2
        add_1 = (2. + mean_red) * diffs[:, 0]
        add_2 = 4. * diffs[:, 1]
        add_3 = (3. + mean_red) * diffs[:, 2]
        return np.sqrt(add_1 + add_2 + add_3)
    get_color_differences = staticmethod(get_color_differences)
