import numpy as np


class FitData:
    """
    Performs a least squares linear fit for given dataset and model to
    determine the source flux(es) and (optionally) blend flux.

    .. _keywords:
    Keywords :
        dataset: :py:class:`~MulensModel.mulensdata.MulensData` object
            A single photometric dataset to be fitted.

        model: :py:class:`~MulensModel.model.Model` object
            The model to fit to the data.

        fix_blend_flux: *False* or *float*
            Default is *False*, i.e. allow the blend flux to be a free
            parameter. If set to a float, it will fix the blend value to that
            value.

        fix_source_flux: *False*, *float*, or *list*
            Default is *False*, i.e. allow the source flux to be a free
            parameter. If set to a float, it will fix the source value to that
            value. For binary source models, a list should be used to set the
            fluxes of the individual sources or fix one and not the other, e.g.
            [2.3, False] would fix f_source_1 to 2.3 but allow a free fit to
            f_source_2. (Not Implemented)

        fix_q_flux: *False* or *float*
            For binary source models, q_flux is the flux ratio between two
            components, i.e. q_flux = f_source_2 / f_source_1
            Default is *False*, i.e. allow the source flux to be a free
            parameter. If set to a float, it will fix the source value to that
            value. (Not Implemented)

    """

    def __init__(self, model=None, dataset=None, fix_blend_flux=False,
                 fix_source_flux=False, fix_q_flux=False):
        self._model = model
        self._dataset = dataset

        # fit parameters
        self.fix_blend_flux = fix_blend_flux
        self.fix_source_flux = fix_source_flux
        self.fix_q_flux = fix_q_flux

        # list containing fluxes of various sources
        self._source_fluxes = None
        self._blend_flux = None

    def _check_for_implementation_errors(self):
        """
        If a setting is not implemented, raise an exception.
        """
        if self.fix_source_flux is not False:
            msg = 'Only fix_source_flux=False is implemented.'
            raise NotImplementedError(msg)

        if self.fix_q_flux is not False:
            msg = 'Only fix_q_flux=False is implemented.'
            raise NotImplementedError(msg)

    def _get_magnifications(self):
        """
        Calculate the model magnifications for the good epochs of the dataset.
        """
        # select = self._dataset.good # Is this faster?

        # currently, model.magnification is good for up to two
        # sources
        if self._model.n_sources == 1:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[self._dataset.good])
        elif self._model.n_sources == 2:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[self._dataset.good], separate=True)
        else:
            msg = ("{0}".format(self._model.n_sources) +
                   " sources used. model.magnification can only" +
                   " handle <=2 sources")
            raise NotImplementedError(msg)

        return mag_matrix

    def _setup_linalg_arrays(self):
        """
        :return: xT and y arrays
        """
        # Find number of fluxes to calculate
        n_fluxes = self._model.n_sources
        n_epochs = np.sum(self._dataset.good)
        y = self._dataset.flux[self._dataset.good]
        x = self._get_magnifications()

        # Account for free or fixed blending
        # Should do a runtime test to compare with lines 83-94
        if self.fix_blend_flux is False:
            x = np.vstack((x, np.ones(n_epochs)))
            n_fluxes += 1
        elif self.fix_blend_flux == 0.:
            pass
        else:
            y -= self.fix_blend_flux

        # Take the transpose of x and define y
        xT = np.copy(x).T
        xT.shape = (n_epochs, n_fluxes)

        # Take into account uncertainties
        sigma_inverse = 1. / self._dataset.err_flux[self._dataset.good]
        y *= sigma_inverse
        xT *= np.array([sigma_inverse] * n_fluxes).T

        return (xT, y)

    def fit_fluxes(self):
        """
        Execute the linear least squares fit to determine the fitted fluxes.
        Sets the values of :py:obj:`~source_fluxes`, :py:obj:`~blend_flux`,
        and (if applicable):py:obj:`~source_flux`.
        """

        self._check_for_implementation_errors()
        (xT, y) = self._setup_linalg_arrays()

        # Solve for the coefficients in y = fs * x + fb (point source)
        # These values are: F_s1, F_s2,..., F_b.
        try:
            results = np.linalg.lstsq(xT, y, rcond=-1)[0]
        except ValueError as e:
            raise ValueError(
                '{0}\n'.format(e) +
                'If either of these numbers ({0}, {1})'.format(
                    np.sum(np.isnan(xT)), np.sum(np.isnan(y))) +
                ' is greater than zero, there is a NaN somewhere,' +
                ' probably in the data.')

        # Record the results
        if self.fix_blend_flux is False:
            self._blend_flux = results[-1]
            self._source_fluxes = results[:-1]
        else:
            self._blend_flux = self.fix_blend_flux
            self._source_fluxes = results

    @property
    def source_flux(self):
        """

        :return: *float*
            the fitted source flux. Only defined for models with a single
            source. See also :py:obj:`~source_fluxes`

            If None, you need to run :py:func:`~fit_fluxes()` to execute the
            linear fit.

        """
        if self._model.n_sources == 1:
            return self.source_fluxes[0]
        else:
            msg = ("source_flux is defined only for models" +
                   " with one source, you have" +
                   " {0}".format(self._model.n_sources) +
                   " sources. Try FitData.source_fluxes instead")

            raise NameError(msg)

    @property
    def source_fluxes(self):
        """

        :return: *np.array*
            the fitted source flux(es).

            If None, you need to run :py:func:`~fit_fluxes()` to execute the
            linear fit.
        """
        return self._source_fluxes

    @property
    def blend_flux(self):
        """

        :return: the fitted blend flux or the value set by
            fix_blend_flux (see :ref:`keywords`).

            If None, you need to run :py:func:`~fit_fluxes()` to execute the
            linear fit.
        """
        return self._blend_flux
