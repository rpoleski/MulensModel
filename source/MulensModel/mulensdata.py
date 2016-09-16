import numpy as np
from astropy.time import Time

class MulensData(object):
    def __init__(self, file_name=None):
        if file_name is not None:
            vector_1, vector_2, vector_3 = np.loadtxt(fname=file_name, unpack=True)
            self._jd_zeropoint = self._get_jd_zeropoint(vector_1)
            self._time = Time(vector_1+self._jd_zeropoint, format="jd")
            self.mag = vector_2
            self.err_mag = vector_3

    @property
    def time(self):
        return self._time.jd - self._jd_zeropoint

    def _get_jd_zeropoint(self, jd_vector):
        """guess what is zeropoint of JD used"""
        if all(jd_vector > 2000.) and all(jd_vector < 10000.):
            return 2450000.
        if all(jd_vector > 52000.) and all(jd_vector < 70000.):
            return 2400000.
        if all(jd_vector > 2452000.):
            return 0.
        raise ValueError('Unrecognized format of JD')

