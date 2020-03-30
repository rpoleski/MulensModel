import matplotlib.pyplot as plt
from MulensModel.utils import Utils

"""Convenience functions for plotting."""


def _get_y_value_y_err(phot_fmt, flux, flux_err):
    """
    just calculate magnitudes if needed, or return input otherwise
    """
    if phot_fmt == 'mag':
        return Utils.get_mag_and_err_from_flux(flux, flux_err)
    else:
        return (flux, flux_err)

def subtract(subtract_2450000=False, subtract_2460000=False):
    """ Check if some value should be subtracted from all dates."""
    subtract = 0.
    if subtract_2450000:
        if subtract_2460000:
            raise ValueError("subtract_2450000 and subtract_2460000 " +
                             "cannot be both True")
        subtract = 2450000.
    if subtract_2460000:
        subtract = 2460000.

    return subtract


def _subtract_xlabel(subtract_2450000, subtract_2460000):
    """
    string that would be past to plt.xlabel()
    """
    # JCY --> plot_functions.py?
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
