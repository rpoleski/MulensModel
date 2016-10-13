import sys
import numpy as np
from astropy.time import Time

from MulensModel.utils import Utils


class MulensData(object):
    def __init__(self, data_list=None, file_name=None, date_fmt="jd", mag_fmt="mag"):
        if data_list is not None and file_name is not None:
            raise ValueError('MulensData cannot be initialized with data_list and file_name')
        elif data_list is not None:
            vector_1, vector_2, vector_3 = list(data_list) # this works for any iterable in input
            self._initialize(date_fmt, mag_fmt, time=vector_1, brightness=vector_2, err_brightness=vector_3)
        elif file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(fname=file_name, unpack=True)
            self._initialize(date_fmt, mag_fmt, time=vector_1, brightness=vector_2, err_brightness=vector_3)
    
    def _initialize(self, date_fmt, mag_fmt, time=None, brightness=None, err_brightness=None):
        """internal function to initialized data using a few numpy arrays"""
        self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt)
        self._time = Time(time+self._date_zeropoint, format="jd")
        self._input_values = brightness
        self._input_values_err = err_brightness        
        self.input_fmt = mag_fmt
        if mag_fmt == "mag":
            self.mag = self._input_values
            self.err_mag = self._input_values_err
            (self.flux, self.err_flux) = Utils.get_flux_and_err_from_mag(mag=self.mag, err_mag=self.err_mag)
        elif mag_fmt == "flux":
            self.flux = self._input_values
            self.err_flux = self._input_values_err
            (self.mag, self.err_mag) = Utils.get_mag_and_err_from_flux(flux=self.flux, err_flux=self.err_flux)
        else:
            raise ValueError('unkonown format of brightness in ' + file_name + ' file')
        self.bad = len(self._time) * [False]

    @property
    def jd(self):
        """full JD time vector"""
        return self._time.jd

    @property
    def time(self):
        """short verion of time vector"""
        return self._time.jd - self._date_zeropoint

    def _get_date_zeropoint(self, date_fmt="jd"):
        """ Return the zeropoint of the date so it can be converted to
        the standard 245#### format."""
        if date_fmt == "jd" or date_fmt == "hjd":
            return 0.
        if date_fmt == "jdprime" or date_fmt == "hjdprime":
            return 2450000.
        if date_fmt == "mjd":
            return 2400000.
        raise ValueError('Invalid value for date_fmt. Allowed values: "jd", "hjd", "jdprime", "hjdprime", "mjd"')

    def _get_jd_zeropoint(self, jd_vector):
        """guess what is zeropoint of JD used"""
        if not hasattr(jd_vector, '__iter__'):
            jd_vector = np.array([jd_vector])
        if all(jd_vector > 2000.) and all(jd_vector < 12000.):
            return 2450000.
        if all(jd_vector > 52000.) and all(jd_vector < 70000.):
            return 2400000.
        if all(jd_vector > 2452000.):
            return 0.
        raise ValueError('Unrecognized format of JD')
