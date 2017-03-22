import numpy as np
from math import fsum

from astropy import __version__ as astropy__version__


MAG_ZEROPOINT = 22. # Defines magnitude at which flux = 1.

month_3letter_to_2digit = {
        'Jan' : '01',
        'Feb' : '02',
        'Mar' : '03',
        'Apr' : '04',
        'May' : '05',
        'Jun' : '06',
        'Jul' : '07',
        'Aug' : '08',
        'Sep' : '09',
        'Oct' : '10',
        'Nov' : '11',
        'Dec' : '12'
        }
        
#JCY: I think this class needs to be a subpackage with related
#functions separated into their own files. e.g. all the flux/mag
#functions together and separated from the math functions.

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
        
    def last_non_space_char_before(line, n_before):
        """find last non-space character before given character"""
        for i in range(n_before-1, 0, -1):
            if line[i] != " ":
                return i
        return -1
    last_non_space_char_before = staticmethod(last_non_space_char_before)

    def complex_fsum(arguments):
        """accurate floating points sum of complex numbers in iterable arguments"""
        real = [arg.real for arg in arguments]
        imag = [arg.imag for arg in arguments]
        return fsum(real) + fsum(imag) * 1j
    complex_fsum = staticmethod(complex_fsum)

    def dot(cartesian, vector):
        """dot product of Astropy CartersianRepresentation and np.array"""
        return cartesian.x * vector[0] + cartesian.y * vector[1] + cartesian.z * vector[2]
    dot = staticmethod(dot)
    
    def combine_dicts(kwargs_to_set, user_dict, default_dict, index=None):
        """construct a new dictionary that that searche user_dict and default_dict
        for keywords from kwargs_to_set list and takes corresponding values,
        if keyword is found
        
        If user_dict[key] is a list or a numpy array and index is provided,
        then user_dict[key][value] is used.
        """
        new_dict = dict()
        for key in kwargs_to_set:
            value = None

            if key in user_dict.keys():
                value = user_dict[key]
            elif key in default_dict.keys():
                value = default_dict[key]
            
            if value is not None:
                if isinstance(value, (list, np.ndarray)) and index is not None:
                    new_dict[key] = value[index]
                else:
                    new_dict[key] = value
        return new_dict
    combine_dicts = staticmethod(combine_dicts)

