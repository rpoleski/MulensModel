"""
Convenience functions for plotting.
"""
from matplotlib.colors import ColorConverter
import numpy as np

from MulensModel.utils import Utils


def _get_y_value_y_err(phot_fmt, flux, flux_err):
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
            'Unrecognized photometry format: {:}, allowed '.format(phot_fmt) +
            'values are "mag" and "flux"')


def subtract(subtract_2450000=False, subtract_2460000=False):
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


def _subtract_xlabel(subtract_2450000=False, subtract_2460000=False):
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
    string that would be past to plt.xlabel()
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


def _color_differences(color_list, color):
    """
    Calculate color difference between a list of colors and a single color.
    Uses algorithm from
    `this Wikipedia page<https://en.wikipedia.org/wiki/Color_difference>`_.

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
