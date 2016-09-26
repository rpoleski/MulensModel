import numpy as np

class Fit(object);
    def __init__(self, data=None, A=None):
        self._datasets = data
        self._A = A #What is A? A-->magnification?
           
    def fit_fluxes(self):
        """fit source(s) and blending fluxes"""
        if len(self._A[0].shape) == 1:
            n_fluxes = 2 #Why have one "fluxes" vector instead of source_flux and blend_flux?
        else:
            n_fluxes, __ = self._A[0].shape 
            n_fluxes += 1
        x = np.empty([n_fluxes, len(dataset)])
        x[n_fluxes-1] = 1.
        for n_data, dataset in enumerate(self._datasets):
            x = np.empty(n_fluxes, len(dataset))
            x[0:(n_fluxes-1),] = self._A[n_data]
            y = data.flux
            results = np.linalg.lstsq(x, y)[0] # F_s1, F_s2..., F_b
        """
        Suggest:
        n_datasets, __ = self._A[0].shape
        for n_data, dataset in enumerate(self._datasets):
            x = np.empty(2.*n_datasets, len(dataset))
            x[0:(2.*n_datasets-1),] = self._A[n_data]
            y = data.flux
            results = np.linalg.lstsq(x, y)[0] # F_s1, F_s2..., F_b
        """
        #How would this work for binary source?

            
class Event(object):
    def __init__(self):
        pass

    def get_chi2(self):
        """calculates chi^2 of current model"""
        self.fit(data=self.datasets, A=self.model.A)
        # self.fit = Fit(data=self.datasets, A=self.model.A) ?
        self.fit.fit_fluxes()
        chi2 = []
        for dataset in self.datasets:
            diff = dataset._input_values - fit.get_input_values(data=dataset)
            select = (not dataset.bad)
            chi2.append((diff[select] / dataset._input_values_err[select])**2)
        self.chi2 = sum(chi2)
        return self.chi2


