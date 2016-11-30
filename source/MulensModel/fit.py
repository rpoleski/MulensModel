import numpy as np
import warnings

from MulensModel.utils import Utils


class Fit(object):
    def __init__(self, data=None, magnification=None, n_sources=None):
        """
        data - a list of MulensModel datasets
        magnification - a list of numpy arrays, 
            each list element is (n_epochs) x (n_sources)
        """
        self._datasets = data 
        self._magnification = magnification
        if magnification is None and n_sources is None:
            raise ValueError('Fit class requires magnifications vectors' 
                    + ' or number of sources directly specified')
        self._n_sources = n_sources
        self._flux_blending = dict()
        self._flux_sources = dict()        
           
    def fit_fluxes(self, fit_blending_all=True):
        """fit source(s) and blending fluxes"""
        if self._n_sources is not None:
            n_sources = self._n_sources
        else:
            if len(self._magnification[0].shape) == 1:
                n_sources = 1
            else:
                n_sources = self._magnification[0].shape[0]
        n_fluxes = n_sources
        if fit_blending_all:
            n_fluxes += 1
        for i_dataset, dataset in enumerate(self._datasets):
            x = np.empty(shape=(n_fluxes, dataset.n_epochs))
            if fit_blending_all:
                x[0:(n_fluxes-1),] = self._magnification[i_dataset]
                x[n_fluxes-1] = 1.
            else:
                x = self._magnification[i_dataset]
            xT = np.copy(x).T
            xT.shape = (dataset.n_epochs, n_fluxes)
            y = np.copy(self._datasets[i_dataset].flux)
            sigma_inverse = 1. / self._datasets[i_dataset].err_flux
            y *= sigma_inverse
            for i, sig_inv in enumerate(sigma_inverse):
                xT[i] *= sig_inv
            
            results = np.linalg.lstsq(xT, y)[0] # These values are: F_s1, F_s2,..., F_b.

            if fit_blending_all:
                self._flux_blending[dataset] = results[-1]
                self._flux_sources[dataset] = results[:-1]
            else:
                self._flux_blending[dataset] = 0.
                self._flux_sources[dataset] = results

    def get_input_format(self, data=None):
        """return model in the same format as given dataset was input"""
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        if self._n_sources is not None:
            n_sources = self._n_sources
        else:
            if len(self._magnification[0].shape) == 1:
                n_sources = 1
            else:
                n_sources = self._magnification[0].shape[0]

        index = self._datasets.index(data)
        if data not in self._flux_blending:
            self._flux_blending[data] = 0. 
            warnings.warn("Blending flux not set. This is strange...", 
                            SyntaxWarning)
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
        
        if data.input_fmt == "mag":
            result = Utils.get_mag_from_flux(flux)
        elif data.input_fmt == "flux":
            result = flux
        else:
            msg = 'Fit.get_input_format() unrecognized data input format'
            raise ValueError(msg)
        return result
    
