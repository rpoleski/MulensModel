import numpy as np


class FitData:
    """
    Fits source and blending fluxes for given data and model magnification.

    Keywords :
        data: :py:class:`MulensData` or *list* of :py:class:`MulensData`
            Photometric data to be fitted.

    """

    def __init__(self, model, dataset, fix_blend_flux=False,
                 fix_source_flux=False, fix_q_flux=False):
        self._model = model
        self._dataset = dataset

        # fit parameters
        self.fix_blend_flux = fix_blend_flux
        self.fix_source_flux = fix_source_flux
        self.fix_q_flux = fix_q_flux

        # list containing fluxes of various sources
        self.source_fluxes = None
        self.blend_flux = None

    def fit_fluxes(self):
        """
        Function to perform
        Arguments:
        Returns itself

        Btw this is too long for a function..
        """

        if self.fix_source_flux is True:
            msg = 'Source flux can only be false'
            raise NotImplementedError(msg)

        n_sources = self._model.n_sources

        # Find number of fluxes to calculate
        n_fluxes = n_sources
        if self.fix_blend_flux is False:
            n_fluxes += 1

        # For dataset, perform a least-squares linear fit for the flux
        # parameters
        dataset = self._dataset
        # suppress bad data
        select = dataset.good
        n_epochs = np.sum(select)

        # Set up the x vector for the linear fit

        x = np.empty(shape=(n_fluxes, n_epochs))

        if self.fix_blend_flux is False:
            # assumes model.magnification returns magnifications for
            # multiple sources

            # currently, model.magnification is good for up to two
            # sources
            if n_sources == 1:
                mag_matrix = self._model.magnification(
                    time=dataset.time[select])
            elif n_sources == 2:
                mag_matrix = self._model.magnification(
                    time=dataset.time[select], separate=True)
            else:
                msg = ("{0}".format(n_sources) +
                       " sources used. model.magnification can only" +
                       " handle <=2 sources")
                raise NotImplementedError(msg)

            # Only deal with good data
            good_mag_matrix = mag_matrix[select]
            x[0:n_sources, ] = good_mag_matrix

            # Row corresponding to blend flux
            x[n_sources] = 1.

        else:
            # fixed blend flux case not implemented yet
            raise NotImplementedError("fixed blend flux case not" +
                                      "implemented")

        # Take the transpose of x and define y
        xT = np.copy(x).T
        xT.shape = (n_epochs, n_fluxes)
        y = dataset.flux[select]
        sigma_inverse = 1. / dataset.err_flux[select]

        y *= sigma_inverse
        xT *= np.array([sigma_inverse] * n_fluxes).T

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
            self.blend_flux = results[-1]
            self.source_fluxes = results[:-1]
        else:
            '''
            self._flux_blending[dataset] = 0.
            self._flux_sources[dataset] = results
            '''
            # fixed blend flux case not implemented yet
            raise NotImplementedError("fixed blend flux case not implemented")

    @property
    def source_flux(self):
        if self._model.n_sources == 1:
            return self.source_fluxes[0]
        else:
            msg = ("source_flux is defined only for models" +
                   " with one source, you have" +
                   " {0}".format(self._model.n_sources) +
                   " sources. Try FitData.source_fluxes")

            raise NameError(msg)
