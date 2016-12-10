import sys
import numpy as np

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.mulenstime import MulensTime

class MulensData(object):
    def __init__(self, data_list=None, file_name=None, date_fmt="jd", 
                 mag_fmt="mag", coords=None):
        date_fmt = date_fmt.lower()
        self._n_epochs = None      
        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            vector_1, vector_2, vector_3 = list(data_list) 
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=coords)
        elif file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(
                fname=file_name, unpack=True, usecols=(0,1,2))
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=coords)
    
    def _initialize(self, date_fmt, mag_fmt, time=None, brightness=None, 
                    err_brightness=None, coords=None):
        """internal function to initialized data using a few numpy arrays"""
        self._time = MulensTime(time=time, date_fmt=date_fmt, coords=coords)
        self._n_epochs = len(time)
        if len(brightness) != self._n_epochs or len(err_brightness) != self._n_epochs:
            raise ValueError('input data in MulesData have different lengths')
        if date_fmt == 'hjd' or date_fmt == 'hjdprime':
            self._time_type = 'hjd'
        else:
            self._time_type = 'jd'
        self._time_corr = None
        self._brightness_input = brightness
        self._brightness_input_err = err_brightness        
        self.input_fmt = mag_fmt
        if mag_fmt == "mag":
            self.mag = self._brightness_input
            self.err_mag = self._brightness_input_err
            (self.flux, self.err_flux) = Utils.get_flux_and_err_from_mag(
                                          mag=self.mag, err_mag=self.err_mag)
        elif mag_fmt == "flux":
            self.flux = self._brightness_input
            self.err_flux = self._brightness_input_err
            (self.mag, self.err_mag) = Utils.get_mag_and_err_from_flux(
                                        flux=self.flux, err_flux=self.err_flux)
        else:
            msg = 'unknown format of brightness in ' + file_name + ' file'
            raise ValueError(msg)
        self.bad = self.n_epochs * [False]

    @property
    def n_epochs(self):
        """give number of epochs"""
        return self._n_epochs

    @property
    def jd(self):
        return self._time.jd

    @property
    def hjd(self):
        """full HJD time vector"""
        return self._time.hjd

    @property
    def time(self):
        """short version of time vector"""
        return self._time.jd - self._time.zeropoint

    @property
    def time_zeropoint(self):
        """return the zeropoint of time vector"""
        return self._time.zeropoint
