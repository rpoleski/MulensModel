import sys
import numpy as np

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u

from MulensModel.utils import Utils
from MulensModel.mulenstime import MulensTime

class MulensData(object):
    def __init__(self, data_list=None, file_name=None, date_fmt="jd", 
                 mag_fmt="mag", coords=None):
        if data_list is not None and file_name is not None:
            m = 'MulensData cannot be initialized with data_list and file_name'
            raise ValueError(m)
        elif data_list is not None:
            vector_1, vector_2, vector_3 = list(data_list) 
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3)
        elif file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(fname=file_name, 
                                                      unpack=True)
            self._initialize(date_fmt, mag_fmt, time=vector_1, 
                             brightness=vector_2, err_brightness=vector_3)

        self._target = None
        if coords is not None:
            if isinstance(coords, SkyCoord):
                self._target = coords
            else:
                msg = 'unsupported format of coords parameter in MulensData()'
                raise ValueError(msg)
    
    def _initialize(self, date_fmt, mag_fmt, time=None, brightness=None, 
                    err_brightness=None):
        """internal function to initialized data using a few numpy arrays"""
        self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt)
        earth_center = EarthLocation.from_geocentric(0., 0., 0., u.m)
        self._time = Time(time+self._date_zeropoint, format="jd", 
                          location=earth_center)
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
        self.bad = len(self._time) * [False]

    @property
    def jd(self):
        """full JD time vector"""
        if self._time_type == 'jd':
            return self._time.jd
        else:
            return (self._time - self._time_correction).jd

    @property
    def hjd(self):
        """full HJD JD time vector"""
        if self._time_type == 'hjd':
            return self._time.jd
        else:
            return (self._time + self._time_correction).jd

    @property
    def _time_correction(self):
        '''time correction: HJD = JD + corr'''
        if self._time_corr is None:
            if self._target is None:
                raise ValueError('Event coordinates in MulensData not set')
            star = SkyCoord(self._target, unit=(u.hour, u.degree), 
                            frame='icrs')
            self._time_corr = self._time.light_travel_time(star, 
                            'heliocentric')
        return self._time_corr

    @property
    def time(self):
        """short verion of time vector"""
        return self._time.jd - self._date_zeropoint

    @property
    def time_zeropoint(self):
        """return the zeropoint of time vector"""
        return self._date_zeropoint

    def _get_date_zeropoint(self, date_fmt="jd"):
        """ Return the zeropoint of the date so it can be converted to
        the standard 245#### format."""
        if date_fmt == "jd" or date_fmt == "hjd":
            return 0.
        if date_fmt == "jdprime" or date_fmt == "hjdprime":
            return 2450000.
        if date_fmt == "mjd":
            return 2400000.5
        lst = '"jd", "hjd", "jdprime", "hjdprime", "mjd"'
        raise ValueError('Invalid value for date_fmt. Allowed values: '+lst)

