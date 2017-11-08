import numpy as np
from math import fsum
from astropy.coordinates import SkyCoord
import astropy.units as u

from MulensModel.utils import Utils
from MulensModel.fit import Fit
from MulensModel.mulensdata import MulensData
from MulensModel.model import Model
from MulensModel.coordinates import Coordinates
            
class Event(object):
    """
    Allows a model to be fit to datasets.

    Arguments :

        :py:obj:`datasets` (required): The data; a
            :py:class:`~MulensModel.mulensdata.MulensData` object or
            list of MulensData objects

        :py:obj:`model` (required): a
            :py:class:`~MulensModel.model.Model` object

        :py:obj:`coords` (optional): the coordinates of the event
            (RA, Dec)

    """
    def __init__(self, datasets=None, model=None, coords=None):
        #Initialise self._model (and check that model is defined).
        if isinstance(model, Model):
            self._model = model
        elif model is None:
            self._model = None
        else:
            raise TypeError('incorrect argument model of class Event()')

        #Initialise self._datasets (and check that datasets is defined).
        if isinstance(datasets, (list, tuple, MulensData)) or datasets is None:
            self._set_datasets(datasets)
        else:
            raise TypeError('incorrect argument datasets of class Event()')

        #Set event coordinates
        if coords is not None:
            self._update_coords(coords=coords)
        else:
            self._coords = None

    @property
    def datasets(self):
        """
        a *list* of :py:class:`~MulensModel.mulensdata.MulensData`
        instances.
        """
        return self._datasets

    @datasets.setter
    def datasets(self, new_value):
        self._set_datasets(new_value)

    @property
    def data_ref(self):
        """
        Reference data set for scaling the model fluxes to (for
        plotting). May be set as a
        :py:class:`~MulensModel.mulensdata.MulensData` object or an
        index (*int*). Default is the first data set.
        """
        return self.model.data_ref
        
    @data_ref.setter
    def data_ref(self, new_value):
        self.model.data_ref = new_value

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

    @property
    def model(self):
        """an instance of :py:class:`~MulensModel.model.Model`"""
        return self._model

    @model.setter
    def model(self, new_value):
        if not isinstance(new_value, Model):
            raise TypeError(('wrong type of Event.model: {:} instead of ' +
                'MulensModel').format(type(new_value)))
        self._model = new_value
        if self._datasets is not None:
            self._model.set_datasets(self._datasets)

        if new_value.coords is not None:
            self._update_coords(coords=new_value.coords)

    @property
    def coords(self):
        """
        see :class:`~MulensModel.coordinates.Coordinates`
        """
        return self._coords
    
    @coords.setter
    def coords(self, new_value):
        self._update_coords(coords=new_value)

    def _update_coords(self, coords=None):
        """Set the coordinates as a SkyCoord object"""
        self._coords = Coordinates(coords)

        if self._model is not None:
            self._model.coords = self._coords

        # We run the command below with try, because _update_coords() is called
        # by _set_datasets before self._datasets is set. 
        try:
            for dataset in self._datasets:
                dataset.coords = self._coords
        except Exception:
            pass

    def get_chi2(self, fit_blending=None):
        """
        Calculates chi^2 of current model by fitting for source and 
        blending fluxes.

        Parameters :
            fit_blending: *boolean*, optional
                If True, then the blend flux is a free parameter. If
                False, the blend flux is fixed at zero.  Default is
                the same as :py:func:`MulensModel.fit.Fit.fit_fluxes()`.

        Returns :
            chi2: *float*
                Chi^2 value

        """
        chi2_per_point = self.get_chi2_per_point(
            fit_blending=fit_blending)
        #Calculate chi^2 given the fit
        chi2 = []
        for i, dataset in enumerate(self.datasets):
            #Calculate chi2 for the dataset excluding bad data 
            select = np.logical_not(dataset.bad)
            chi2.append(fsum(chi2_per_point[i][select]))

        self.chi2 = fsum(chi2)
        return self.chi2

    def get_chi2_per_point(self, fit_blending=None):
        """Calculates chi^2 for each data point of the current model by
        fitting for source and blending fluxes.

        Parameters :
            fit_blending: *boolean*, optional
                Are we fitting for blending flux? If not then blending flux is 
                fixed to 0.  Default is the same as
                :py:func:`MulensModel.fit.Fit.fit_fluxes()`.

        Returns :
            chi2: *np.ndarray*  
                Chi^2 contribution from each data point,
                e.g. chi2[obs_num][k] returns the chi2 contribution
                from the *k*-th point of observatory *obs_num*.
        """
        #Define a Fit given the model and perform linear fit for fs and fb
        self.fit = Fit(data=self.datasets, 
                       magnification=self.model.data_magnification) 
        if fit_blending is not None:
            self.fit.fit_fluxes(fit_blending=fit_blending)
        else:
            self.fit.fit_fluxes()

        #Calculate chi^2 given the fit
        chi2_per_point = []
        for i, dataset in enumerate(self.datasets):
            diff = dataset._brightness_input \
                 - self.fit.get_input_format(data=dataset)
            chi2_per_point.append(
                (diff / dataset._brightness_input_err)**2)

        chi2_per_point = np.array(chi2_per_point)
        return chi2_per_point


    def get_ref_fluxes(self, data_ref=None):
        """
        Get source and blending fluxes for the reference dataset. See 
        :py:func:`MulensModel.model.Model.get_ref_fluxes()` for details.
        """
        return self.model.get_ref_fluxes(data_ref=data_ref)

    def plot_model(self, **kwargs):
        """
        Plot the model light curve in magnitudes. See
        :py:func:`MulensModel.model.Model.plot_lc()` for details.
        """
        self.model.plot_lc(**kwargs)

    def plot_data(self, **kwargs):
        """
        Plot the data scaled to the model. See 
        :py:func:`MulensModel.model.Model.plot_data()` for details.
        """
        self.model.plot_data(**kwargs)

    def plot_residuals(self, **kwargs):
        """
        Plot the residuals (in magnitudes) of the model. 
        See :py:func:`MulensModel.model.Model.plot_residuals()` for details.
        """
        self.model.plot_residuals(**kwargs)

    def clean_data(self):
        """masks outlying datapoints. **Not Implemented.**"""
        raise NotImplementedError("This feature has not been implemented yet")

    def estimate_model_params(self):
        """estiamtes model parameters without fitting them. 
        **Not Implemented.**"""
        raise NotImplementedError("This feature has not been implemented yet")

