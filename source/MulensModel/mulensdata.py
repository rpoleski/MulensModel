import sys
import numpy as np
from astropy.time import Time

class MulensData(object):
    def __init__(self, file_name=None,date_fmt="jd"):
        if file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(fname=file_name, unpack=True)
            self._date_zeropoint = self._get_date_zeropoint(date_fmt=date_fmt)
            self._time = Time(vector_1+self._date_zeropoint, format="jd")
            self.mag = vector_2
            self.err_mag = vector_3

    @property
    def jd(self):
        return self._time.jd

    @property
    def time(self):
        return self._time.jd - self._date_zeropoint

    def _get_date_zeropoint(self,date_fmt="jd"):
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

