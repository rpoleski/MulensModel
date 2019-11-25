import numpy as np
import warnings

from MulensModel.utils import Utils


class Fit(object):
    """
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

    def __init__(self, data=None, magnification=None, n_sources=None):
        # Initialize self._datasets, self._magnification, and self._n_sources
        if isinstance(data, list):
            self._datasets = data
        else:
            self._datasets = [data]

        if isinstance(magnification, list):
            self._magnification = magnification
        else:
            self._magnification = [magnification]

        if magnification is None and n_sources is None:
            raise ValueError(
                'Fit class requires magnifications vectors' +
                ' or number of sources directly specified')
        self._n_sources = n_sources

        # Set up numpy ndarrays for flux parameters
        self._flux_blending = dict()
        self._flux_sources = dict()

    def fit_fluxes(self, fit_blending=True):
        """
        Fit source(s) and blending fluxes. Performs a least squares
        linear fit (*np.linalg.lstsq()*) to the data for the flux
        parameters. I.e., given the data :math:`y` and magnification
        :math:`A`, solves for :math:`f_{source}` and
        :math:`f_{blend}`:

        .. math::

           y = f_{source} * A + f_{blend}

        Parameters :
            fit_blending: *boolean*, optional
                Do you want the blending to be fitted (*True*) or to
                be fixed at 0 (*False*)? Default is *True*

        """
        if fit_blending not in {True, False}:
            msg = 'Unrecognized format of fit_blending: {:}'
            raise ValueError(msg.format(fit_blending))
        n_sources = self.get_n_sources()

        # Add parameters for blended light (if appropriate)
        n_fluxes = n_sources
        if fit_blending:
            n_fluxes += 1

        # For each dataset, perform a least-squares linear fit for the flux
        # parameters
        for (i_dataset, dataset) in enumerate(self._datasets):
            # suppress bad data
            select = dataset.good
            n_epochs = np.sum(select)
            all_ok = np.all(select)

            # Set up the x vector for the linear fit
            x = np.empty(shape=(n_fluxes, n_epochs))
            if fit_blending:
                if all_ok:
                    x[0:n_sources, ] = self._magnification[i_dataset]
                else:
                    x[0:n_sources, ] = self._magnification[i_dataset][select]
                x[n_sources] = 1.
            else:
                if all_ok:
                    x = self._magnification[i_dataset]
                else:
                    x = self._magnification[i_dataset][select]

            # Take the transpose of x and define y
            xT = np.copy(x).T
            xT.shape = (n_epochs, n_fluxes)
            if all_ok:
                y = np.copy(self._datasets[i_dataset].flux)
                sigma_inverse = 1. / self._datasets[i_dataset].err_flux
            else:
                y = self._datasets[i_dataset].flux[select]
                sigma_inverse = 1. / self._datasets[i_dataset].err_flux[select]
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
                    ' probably in the data (dataset={0}).'.format(i_dataset))

            # Record the results
            if fit_blending:
                self._flux_blending[dataset] = results[-1]
                self._flux_sources[dataset] = results[:-1]
            else:
                self._flux_blending[dataset] = 0.
                self._flux_sources[dataset] = results

    def flux_of_sources(self, dataset):
        """
        Fluxes of source(s).

        Parameters :
            dataset: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which source fluxes will be given.

        Returns :
            source_fluxes: *np.ndarray*
                Fluxes of sources in units where 1 corresponds to 22 mag.
                The number of array elements is the same as number of sources.

        """
        return self._flux_sources[dataset]

    def blending_flux(self, dataset):
        """
        Blending flux for given dataset.

        Parameters :
            dataset: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which blending flux will be given.

        Returns :
            blending_flux: *np.float64*
                Blending flux in units where 1 corresponds to 22 mag.

        """
        return self._flux_blending[dataset]

    def get_flux(self, data):
        """
        Microlensing model in flux units.

        Parameters :
            data: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which model will be returned.

        Returns :
            flux: *np.ndarray*
                Microlensing model in flux units.
        """
        if data not in self._flux_blending:
            self._flux_blending[data] = 0.
            warnings.warn(
                "Blending flux not set. This is strange...", SyntaxWarning)

        n_sources = self.get_n_sources()
        index = self._datasets.index(data)

        flux = np.ones(data.n_epochs) * self._flux_blending[data]
        if n_sources == 1:
            if data not in self._flux_sources:
                self._flux_sources[data] = 1.
                warnings.warn(
                    "Source flux not set. This is strange...", SyntaxWarning)
            flux += self._flux_sources[data] * self._magnification[index]
        else:
            for i in range(n_sources):
                flux += self._flux_sources[data][i] \
                      * self._magnification[index][i]

        return flux

    def get_input_format(self, data=None):
        """
        Microlensing model in the same format as given dataset. The
        output is either in flux units or magnitudes, depending on
        format of the input data.

        Parameters :
            data: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which model will be returned.

        Returns :
            model: *np.ndarray*
                Microlensing model in flux units or magnitudes (depending on
                the format of input data).

        """
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        if data not in self._flux_blending:
            self._flux_blending[data] = 0.
            warnings.warn(
                "Blending flux not set. This is strange...", SyntaxWarning)

        # Return the model flux in either flux or magnitudes
        # (depending on data)
        if data.input_fmt == "mag":
            result = Utils.get_mag_from_flux(self.get_flux(data))
        elif data.input_fmt == "flux":
            result = self.get_flux(data)
        else:
            msg = 'Fit.get_input_format() unrecognized data input format'
            raise ValueError(msg)
        return result

    def get_chi2_format(self, data):
        """
        Microlensing model in the format used for chi^2 calculation.
        The output is in flux space in most cases, but can be in
        magnitudes depending on dataset format.

        Parameters :
            data: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which model will be returned.

        Returns :
            model: *np.ndarray*
                Microlensing model in flux units or magnitudes (depending on
                the settings of input data).

        """
        if data not in self._flux_blending:
            self._flux_blending[data] = 0.
            warnings.warn(
                "Blending flux not set. This is strange...", SyntaxWarning)

        if data.chi2_fmt == "mag":
            result = Utils.get_mag_from_flux(self.get_flux(data))
        elif data.chi2_fmt == "flux":
            result = self.get_flux(data)
        else:
            msg = 'Fit.get_chi2_format() unrecognized data input format'
            raise ValueError(msg)
        return result

    def get_n_sources(self):
        """
        Count sources.

        Returns :
            n_sources: *int*
                Number of sources in input data.
        """

        # if not set, determine from the shape of _magnification
        if self._n_sources is not None:
            n_sources = self._n_sources
        else:
            if len(self._magnification[0].shape) == 1:
                n_sources = 1
            else:
                n_sources = self._magnification[0].shape[0]
        return n_sources

    def update(self, fit):
        """
        Update internal variables using information from a different
        instance of the same class.

        Parameters :
            fit: :py:class:`~MulensModel.fit.Fit`
                A different instance of this class.
        """
        if not isinstance(fit, Fit):
            raise TypeError('fit parameter must be of MulensModel.Fit type')
        if self.get_n_sources() != fit.get_n_sources():
            raise ValueError('internal erro with number of sources')

        self._datasets.extend(fit._datasets)
        self._magnification.extend(fit._magnification)
        self._flux_blending.update(fit._flux_blending)
        self._flux_sources.update(fit._flux_sources)
