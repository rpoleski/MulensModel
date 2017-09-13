import numpy as np
import warnings

from MulensModel.utils import Utils


class Fit(object):
    """
    fit_fluxes(): Performs a least squares linear fit
        (np.linalg.lstsq) to the data for the flux
        parameters. i.e. given the data y and magnifcation A, solves
        for f_source and f_blend:
    
        y = f_source * A + f_blend

        Allows for zero blending.

    get_input_format(data): Returns the model magnifications in the
        same flux system as data.
    """

    def __init__(self, data=None, magnification=None, n_sources=None):
        """
        Args:
            data: a list of MulensModel.MulensData objects
            magnification: the model magnification for each data epoch: a list
                of numpy arrays, each list element is (n_epochs) x (n_sources)
            n_sources (optional): the number of microlens sources
        """
        # Initialize self._datasets, self._magnification, and self._n_sources
        if isinstance(data, list):
            self._datasets = data 
        else:
            self._datasets = [data]

        if isinstance(magnification, list):
            self._magnification = magnification
        else:
            self._magnificaiton = [magnification]

        if magnification is None and n_sources is None:
            raise ValueError('Fit class requires magnifications vectors' 
                    + ' or number of sources directly specified')
        self._n_sources = n_sources
        
        #Set up numpy arrays for flux parameters
        self._flux_blending = dict()
        self._flux_sources = dict()

           
    def fit_fluxes(self, fit_blending=True):
        """fit source(s) and blending fluxes"""
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
                x[0:(n_fluxes-1),] = self._magnification[i_dataset][select]
                x[n_fluxes-1] = 1.
            else:
                x = self._magnification[i_dataset][select]

            # Take the transpose of x and define y
            xT = np.copy(x).T
            xT.shape = (n_epochs, n_fluxes)
            y = np.copy(self._datasets[i_dataset].flux[select])
            sigma_inverse = 1. / self._datasets[i_dataset].err_flux[select]
            y *= sigma_inverse
            for (i, sig_inv) in enumerate(sigma_inverse):
                xT[i] *= sig_inv

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
        """return blending flux for given dataset"""
        return self._flux_blending[dataset]
        
    def flux_of_sources(self, dataset):
        """return fluxes for all sources in a numpy array"""
        return self._flux_sources[dataset]

    def get_input_format(self, data=None):
        """
        return model in the same format as data, i.e. use data as the
        reference dataset for the flux system. data must be in the
        list of datasets.
        """
        #Check for potential problems
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        if data not in self._flux_blending:
            self._flux_blending[data] = 0. 
            warnings.warn("Blending flux not set. This is strange...", 
                            SyntaxWarning)

        #Initialize parameters
        n_sources = self.get_n_sources()
        index = self._datasets.index(data)

        #Calculate the model flux
        flux = np.ones(data.n_epochs) * self._flux_blending[data]
        if n_sources == 1:
            if data not in self._flux_sources:
                self._flux_sources[data] = 1. 
                warnings.warn("Source flux not set. This is strange...", 
                                SyntaxWarning)
            flux += self._flux_sources[data] * self._magnification[index]
        else:
            for i in range(n_sources):
                flux += self._flux_sources[data][i] \
                      * self._magnification[index][i]
        
        #Return the model flux in either flux or magnitudes (depending on data)
        if data.input_fmt == "mag":
            result = Utils.get_mag_from_flux(flux)
        elif data.input_fmt == "flux":
            result = flux
        else:
            msg = 'Fit.get_input_format() unrecognized data input format'
            raise ValueError(msg)
        return result

    
    def get_n_sources(self):
        """count n_sources (number of sources)"""
        #if not set, determine from the shape of _magnification
        if self._n_sources is not None:
            n_sources = self._n_sources
        else:
            if len(self._magnification[0].shape) == 1:
                n_sources = 1
            else:
                n_sources = self._magnification[0].shape[0]
        return n_sources
