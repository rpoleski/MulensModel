import sys
import numpy as np

from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.mulenstime import MulensTime

class MulensData(object):
    def __init__(self, data_list=None, file_name=None, date_fmt="jd", 
                 mag_fmt="mag", coords=None, ra=None, dec=None):
        date_fmt = date_fmt.lower()
        self._n_epochs = None   

        coords_msg = 'Must specify both or neither of ra and dec'
        self._coords = None
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._coords = coords
            else:
                self._coords = SkyCoord(coords, unit=(u.hourangle, u.deg))
        if ra is not None:
            if dec is not None:
                self._coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            else:
                raise AttributeError(coords_msg)
        else:
            if ra is not None:
                raise AttributeError(coords_msg)

        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            vector_1, vector_2, vector_3 = list(data_list) 
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
        elif file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(
                fname=file_name, unpack=True, usecols=(0,1,2))
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3,
                             coords=self._coords)
    
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

    @property
    def coords(self):
        """
        Sky coordinates (RA,Dec)
        """
        return self._coords

    @coords.setter
    def coords(self, new_value):
        if isinstance(new_value, SkyCoord):
            self._coords = new_value
        else:
            self._coords = SkyCoord(new_value, unit=(u.hourangle, u.deg))
        self._time._target = self._coords

    @property
    def ra(self):
        """
        Right Ascension
        """
        return self._coords.ra

    @ra.setter
    def ra(self, new_value):
        try:
            self._coords.ra = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    new_value, 0.0, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    new_value, self._coords.dec, unit=(u.hourangle, u.deg)) 
        self._time._target = self._coords

    @property
    def dec(self):
        """
        Declination
        """
        return self._coords.dec

    @dec.setter
    def dec(self, new_value):
        try:
            self._coords.dec = new_value
        except AttributeError:
            if self._coords is None:
                self._coords = SkyCoord(
                    0.0, new_value, unit=(u.hourangle, u.deg))
            else:
                self._coords = SkyCoord(
                    self._coords.ra, new_value, unit=(u.hourangle, u.deg))
        self._time._target = self._coords
