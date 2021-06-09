import numpy as np
import warnings

from MulensModel.utils import Utils


class Fit(object):
    """
    DEPRECATED

    Fits source and blending fluxes for given data and model magnification.

    Keywords :
        data: :py:class:`MulensData` or *list* of :py:class:`MulensData`
            Photometric data to be fitted.

        magnification: *np.ndarray* or *list of np.ndarrays*
            Model magnification.

        n_sources: *int*
            The number of microlensing sources. *It's suggested not to use
            this option now.*

    """

    def __init__(self, **kwargs):
        raise NameError(
            'The Fit() class has been deprecated. The overall architecture ' +
            'of MulensModel has been reworked. Hopefully, you never used ' +
            'the Fit() class directly, but if you did, please see Event() ' +
            'and FitData(). ')
