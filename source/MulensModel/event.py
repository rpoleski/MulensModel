import numpy as np
from math import fsum
from astropy.coordinates import SkyCoord
import astropy.units as u

from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData
from MulensModel.model import Model

            
class Event(object):
    """
    Connects datasets to a model.

    Attributes:
        datasets: the input data (MulensModel.MulensData)
        model: the microlensing model (MulensModel.Model)

        coords: the sky coordinates
        ra: Right Ascension
        dec: Declination
    """
    def __init__(self, datasets=None, model=None, coords=None):
        """
        Create an Event object.

        Args:
            datasets (required): The data; a MulensData object or list of 
                MulensData objects
            model (required): a MulensModel.Model object
            coords (optional): the coordinates of the event        
        """
        #Initialize self._model (and check that model is defined)
        if isinstance(model, Model):
            self._model = model
        elif model is None:
            self._model = None
        else:
            raise TypeError('incorrect argument model of class Event()')

        #Initialize self._datasets (and check that datasets is defined)
        if isinstance(datasets, (list, tuple, MulensData)) or datasets is None:
            self._set_datasets(datasets)
        else:
            raise TypeError('incorrect argument datasets of class Event()')

        #Set event coordinates
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
        """
        sets the value of self._datasets
        can be called by __init__ or @datasets.setter
        passes datasets to property self._model
        """
        if isinstance(new_value, list):
            for dataset in new_value:
                if dataset.coords is not None:
                    self._update_coords(coords=dataset.coords)
        if isinstance(new_value, MulensData):
            if new_value.coords is not None:
                self._update_coords(coords=new_value.coords)
            new_value = [new_value]
        if new_value is None:
            self._datasets = None
            return
        self._datasets = new_value
        if isinstance(self._model, Model):
            self._model.set_datasets(self._datasets)


    def get_chi2(self, fit_blending_all=None):
        """calculates chi^2 of current model by fitting for fs and fb"""
        #Define a Fit given the model and perform linear fit for fs and fb
        self.fit = Fit(data=self.datasets, 
                       magnification=self.model.data_magnification) 
        if fit_blending_all is not None:
            self.fit.fit_fluxes(fit_blending_all=fit_blending_all)
        else:
            self.fit.fit_fluxes()

        #Calculate chi^2 given the fit
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
        """Set the coordinates as a SkyCoord object"""
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

    def plot_model(self, 
        times=None, t_range=None, t_start=None, t_stop=None, dt=None, 
        n_epochs=None, data_ref=None, f_source=None, f_blend=None, **kwargs):
        """
        Plot the model lightcurve in magnitudes scaled to data_ref
        (either an index or a MulensData object). If data_ref is not
        specified or data_ref is None, it will use the first dataset
        (see Model.get_ref_fluxes).
        """
        self.model.plot_lc( 
            times=times, t_range=t_range, t_start=t_start, t_stop=t_stop, 
            dt=dt, n_epochs=n_epochs, data_ref=data_ref, f_source=f_source, 
            f_blend=f_blend,**kwargs)

    def plot_data(self, data_ref=None, errors=True, **kwargs):
        """
        Plot the data scaled to the same flux system specified by
        data_ref. Uses the model to calculate the magnifications.
        """
        self.model.plot_data(data_ref=data_ref, errors=errors, **kwargs)

    def plot_residuals(self, errors=True, **kwargs):
        self.model.plot_residuals(errors=errors, **kwargs)
    

if __name__ == "__main__":
    #Generate a model
    t_0 = 5380.
    u_0 = 0.5
    t_E = 18.
    model = Model(t_0=t_0, u_0=u_0, t_E=t_E)
    
    #Generate fake data offset from that model
    times = np.arange(5320,5420.)

    f_source = 0.1
    f_blend = 0.5
    mod_fluxes = f_source * model.magnification(times) + f_blend

    I_mag = Utils.get_mag_from_flux(mod_fluxes)

    random_numbers = np.random.randn(times.size)
    errors = 0.01 * np.ones(times.shape)
    
    mags = I_mag + random_numbers * errors

    data = MulensData(data_list=[times, mags, errors])

    #Generate event and fit
    event = Event(datasets=data, model=model)

    print(event.get_chi2(), np.sum(np.abs(random_numbers)**2))
    
    import matplotlib.pyplot as pl
    pl.errorbar(times, mags, yerr=errors)
    model.plot_lc(times=times, data_ref=data)
    event.model.plot_lc(times=times, data_ref=data, linestyle='--')

    pl.show()
