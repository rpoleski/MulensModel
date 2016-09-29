import numpy as np

from MulensModel.utils import Utils


class Fit(object);
    def __init__(self, data=None, magnification=None):
        self._datasets = data
        self._magnification = magnification
           
    def fit_fluxes(self):
        """fit source(s) and blending fluxes"""
        if len(self._magnification[0].shape) == 1:
            n_sources = 1
        else:
            n_sources = self._magnification[0].shape[0]
        n_fluxes = n_sources + 1
        self._flux_blending = dict()
        self._flux_sources = dict()
        for i_dataset, dataset in enumerate(self._datasets):
            x = np.empty(n_fluxes, len(dataset))
            x[0:(n_fluxes-1),] = self._magnification[i_dataset]
            x[n_fluxes-1] = 1.
            results = np.linalg.lstsq(x, self._datasets[i_dataset])[0] # F_s1, F_s2..., F_b
            self._flux_blending[dataset] = results[-1]
            self._flux_sources[dataset] = result[:-1]

    def get_input_format(self, data=None):
        """return model in the same format as given dataset was input"""
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        n_sources = self._magnification[0].shape[0]
        flux = np.ones() * self._flux_blending[data]

        if len(self._magnification[0].shape) == 1:
            n_sources = 1
        else:
            n_sources = self._magnification[0].shape[0]
        index = self._datasets.index(data)

        if n_sources == 1:
            flux += self._flux_sources[data][0] * self._magnification[index]
        else
            for i in range(n_sources):
                flux += self._flux_sources[data][i] * self._magnification[index][i]

        if data.mag_fmt == "mag":
            result = get_mag_from_flux(flux)
        elif data.mag_fmt == "flux":
            result = flux
        else:
            raise ValueError('Fit.get_input_format() unrecognized data input format')
        return result
    
            

