import numpy as np

from MulensModel.utils import Utils


class Fit(object):
    def __init__(self, data=None, magnification=None):
        """
        data - a list of MulensModel datasets
        magnification - a list of numpy arrays, each list element is (n_epochs) x (n_sources)
        """
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
            x = np.empty(shape=(n_fluxes, len(dataset.jd)))
            x[0:(n_fluxes-1),] = self._magnification[i_dataset]
            x[n_fluxes-1] = 1.
            xT = x.T
            y = self._datasets[i_dataset].flux
            variance_inverse = self._datasets[i_dataset].err_flux**-2
            y *= variance_inverse
            for i,vari in enumerate(variance_inverse):
                xT[i] *= vari
            results = np.linalg.lstsq(xT, y)[0] # F_s1, F_s2..., F_b
            #results[0] = 1664.6791620618908
            #results[1] = 27.118662350933683

            # print(sum(variance_inverse * (np.dot(xx.T, results) - self._datasets[i_dataset].flux)**2))
            #print(results) # 41.81485 0.68119 for MAG_ZEROPOINT = 18
            self._flux_blending[dataset] = results[-1]
            self._flux_sources[dataset] = results[:-1]

    def get_input_format(self, data=None):
        """return model in the same format as given dataset was input"""
        if data is None:
            raise ValueError('Fit.get_input_format() dataset not provided')
        if len(self._magnification[0].shape) == 1:
            n_sources = 1
        else:
            n_sources = self._magnification[0].shape[0]            

        index = self._datasets.index(data)

        flux = np.ones(len(self._magnification[0])) * self._flux_blending[data]

        if n_sources == 1:
            flux += self._flux_sources[data] * self._magnification[index]
        else:
            for i in range(n_sources):
                flux += self._flux_sources[data][i] * self._magnification[index][i]

        if data.input_fmt == "mag":
            result = Utils.get_mag_from_flux(flux)
        elif data.input_fmt == "flux":
            result = flux
        else:
            raise ValueError('Fit.get_input_format() unrecognized data input format')
        return result
    
