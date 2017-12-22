import numpy as np
import warnings

from MulensModel.utils import Utils


class Fit(object):
    """
    Fits source and blending fluxes for given data and model magnification.

    Keywords :
        data: :py:class:`MulensData` or list of :py:class:`MulensData`
        instances
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
        n_sources = self.get_n_sources()

        # Add parameters for blended light (if appropriate)
        n_fluxes = n_sources
        if fit_blending:
            n_fluxes += 1

        # For each dataset, perform a least-squares linear fit for the flux
        # parameters
        for (i_dataset, dataset) in enumerate(self._datasets):
            # suppress bad data
            select = np.logical_not(dataset.bad)
            n_epochs = dataset.n_epochs - np.sum(dataset.bad)

            # Set up the x vector for the linear fit
            x = np.empty(shape=(n_fluxes, n_epochs))
            if fit_blending:
                x[0:(n_fluxes-1), ] = self._magnification[i_dataset][select]
                x[n_fluxes-1] = 1.
            else:
                x = self._magnification[i_dataset][select]

            # Take the transpose of x and define y
            xT = np.copy(x).T
            xT.shape = (n_epochs, n_fluxes)
            y = np.copy(self._datasets[i_dataset].flux[select])
            sigma_inverse = 1. / self._datasets[i_dataset].err_flux[select]
            y *= sigma_inverse
            if fit_blending:
                xT *= np.array([sigma_inverse, sigma_inverse]).T
            else:
                xT *= np.array([sigma_inverse]).T

            # Solve for the coefficients in y = fs * x + fb (point source)
            # These values are: F_s1, F_s2,..., F_b.
            results = np.linalg.lstsq(xT, y)[0]

            # Record the results
            if fit_blending:
                self._flux_blending[dataset] = results[-1]
                self._flux_sources[dataset] = results[:-1]
            else:
                self._flux_blending[dataset] = 0.
                self._flux_sources[dataset] = results

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

    def get_input_format(self, data=None):
        """
        Microlensing model in the same format as given dataset. The
        output is either in flux units or magnitudes, depending on
        format of the input data.

        Parameters :
            data: :py:class:`~MulensModel.mulensdata.MulensData`
                A dataset for which model will be returned.

        Returns :
            magnification: *np.ndarray*
                Magnifications in flux units or magnitudes (depending on
                the format of input data).

        """
        # Check for potential problems
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        if data not in self._flux_blending:
            self._flux_blending[data] = 0.
            warnings.warn(
                "Blending flux not set. This is strange...", SyntaxWarning)

        # Initialize parameters
        n_sources = self.get_n_sources()
        index = self._datasets.index(data)

        # Calculate the model flux
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

        # Return the model flux in either flux or magnitudes
        # (depending on data)
        if data.input_fmt == "mag":
            result = Utils.get_mag_from_flux(flux)
        elif data.input_fmt == "flux":
            result = flux
        else:
            msg = 'Fit.get_input_format() unrecognized data input format'
            raise ValueError(msg)
        return result

    def get_n_sources(self):
        """Count sources.

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
