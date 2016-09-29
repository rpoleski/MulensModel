import numpy as np

from MulensModel.utils import Utils
from MulensModel.fit import Fit

            
class Event(object):
    def __init__(self):
        pass

    def get_chi2(self):
        """calculates chi^2 of current model"""
        self.fit = Fit(data=self.datasets, magnification=self.model.magnification) 
        self.fit.fit_fluxes()
        chi2 = []
        for dataset in self.datasets:
            diff = dataset._input_values - fit.get_input_format(data=dataset)
            select = (not dataset.bad)
            chi2.append((diff[select] / dataset._input_values_err[select])**2)
        self.chi2 = sum(chi2)
        return self.chi2


