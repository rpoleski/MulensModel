import numpy as np
from math import fsum
from astropy.coordinates import SkyCoord
import astropy.units as u

from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData
from MulensModel.model import Model

            
class Event(object):
    def __init__(self, datasets=None, model=None, coords=None):
        if isinstance(model, Model):
            self._model = model
        elif model is None:
            self._model = None
        else:
            raise TypeError('incorrect argument model of class Event()')
        if isinstance(datasets, (list, tuple, MulensData)) or datasets is None:
            self._set_datasets(datasets)
        else:
            raise TypeError('incorrect argument datasets of class Event()')
        if coords is not None:
            self._update_coords(coords=coords)

    @property
    def datasets(self):
        """a list of MulensData instances that represent all event datasets"""
        return self._datasets

    @datasets.setter
    def datasets(self, new_value):
        self._set_datasets(new_value)

    @property
    def model(self):
        """an instance of Model"""
        return self._model

    @model.setter
    def model(self, new_value):
        #Needs a check for MulensModel class
        self._model = (new_value)
        try:
            if self._datasets is not None:
                self._model.set_datasets(self._datasets)
        except:
            pass
        if new_value.coords is not None:
            self._update_coords(coords=new_value.coords)

    def _set_datasets(self, new_value):
        """sets the value of self._datasets
        can be called by __init__ or @datasets.setter
        passes datasets to property self._model"""
        if isinstance(new_value, list):
            for dataset in new_value:
                if dataset.coords is not None:
                    self._update_coords(coords=dataset.coords)
        if isinstance(new_value, MulensData):
            self._update_coords(coords=new_value.coords)
            new_value = [new_value]
        if new_value is None:
            self._datasets = None
            return
        self._datasets = new_value
        if isinstance(self._model, Model):
            self._model.set_datasets(self._datasets)

    def get_chi2(self, fit_blending_all=None):
        """calculates chi^2 of current model"""
        self.fit = Fit(data=self.datasets, 
                       magnification=self.model.magnification) 
        if fit_blending_all is not None:
            self.fit.fit_fluxes(fit_blending_all=fit_blending_all)
        else:
            self.fit.fit_fluxes()
        chi2 = []
        for dataset in self.datasets:
            diff = dataset._brightness_input \
                 - self.fit.get_input_format(data=dataset)
            select = np.logical_not(dataset.bad)
            chi2.append(fsum((diff[select] 
                        / dataset._brightness_input_err[select])**2))
        self.chi2 = fsum(chi2)
        return self.chi2

    def clean_data(self):
        """masks outlying datapoints"""
        pass

    def estimate_model_params(self):
        """estiamtes model parameters without fitting them"""
        pass

    @property
    def coords(self):
        "Return the event Sky Coordinates (RA, Dec)"
        return self._coords
    
    @coords.setter
    def coords(self, new_value):
        self._update_coords(coords=new_value)

    @property
    def ra(self):
        """
        Right Ascension
        """
        return self._coords.ra

    @ra.setter
    def ra(self, new_value):
        pass
        """
        try:
            self._coords.ra = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    new_value, 0.0, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    new_value, self._coords.dec, unit=(u.hourangle, u.deg)) 
        """

    @property
    def dec(self):
        """
        Declination
        """
        return self._coords.dec

    @dec.setter
    def dec(self, new_value):
        pass
        """
        try:
            self._coords.dec = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    0.0, new_value, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    self._coords.ra, new_value, unit=(u.hourangle, u.deg))
        """


    def _update_coords(self, coords=None):
        if isinstance(coords, SkyCoord):
            self._coords = coords
        else:
            self._coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        try:
            self._model.coords = self._coords
        except:
            pass
        try:
            for dataset in self._datasets:
                dataset.coords = self._coords
        except:
            pass
    
