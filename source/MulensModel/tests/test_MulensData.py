import sys
import unittest
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

from MulensModel.mulensdata import MulensData


for path in sys.path:
    if path.find("MulensModel") > 0:
        MODULE_PATH = "/".join(path.split("/")[:-1])
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"


def test_file_read():
    '''read sample file and check if values match'''
    data = MulensData(file_name=SAMPLE_FILE_01)

    np.testing.assert_almost_equal(data.time[0], 5264.84100, 
                                   err_msg="time of first line doesn't match")
    
    assert data.mag[-1] == 13.913, "magnitude of the last line doesn't match"

def test_HJD_JD_conversion():
    '''test if heliocentric correction is calculated properly'''
    t_jd = np.array([7572.458333])
    t_hjd = np.array([7572.464055])
    mag = 16. + 0. * t_jd
    err = 0.01 + 0. * t_jd
    coords = SkyCoord('18:00:00 -30:00:00', unit=(u.hourangle, u.deg))
    d_jd = MulensData([t_jd, mag, err],  date_fmt="jdprime", coords=coords)
    np.testing.assert_almost_equal(d_jd.hjd-d_jd.time_zeropoint, t_hjd)
    d_hjd =  MulensData([t_hjd, mag, err],  date_fmt="hjdprime", coords=coords)
    np.testing.assert_almost_equal(d_jd.jd-d_jd.time_zeropoint, t_jd)

def test_get_date_zeropoint_1():
    test_data = MulensData()
    assert test_data._get_date_zeropoint(date_fmt="jd") == 0.

def test_get_date_zeropoint_2():
    test_data = MulensData()
    assert test_data._get_date_zeropoint(date_fmt="hjd") == 0.

def test_get_date_zeropoint_3():
    test_data = MulensData()
    assert test_data._get_date_zeropoint(date_fmt="jdprime") == 2450000.

def test_get_date_zeropoint_4():
    test_data = MulensData()
    assert test_data._get_date_zeropoint(date_fmt="hjdprime") == 2450000.

def test_get_date_zeropoint_5():
    test_data = MulensData()
    assert test_data._get_date_zeropoint(date_fmt="mjd") == 2400000.5
    
def test_data_list_1():
    t = np.array([7500., 7501.])
    m = np.array([21.0, 21.1])
    e = np.array([0.001, 1.000])
    data = MulensData([t, m, e], date_fmt="jdprime")
    np.testing.assert_almost_equal(data.time, t, err_msg
                                   ='problem with time vector in MulensData')
    np.testing.assert_almost_equal(data.time_zeropoint, 2450000., 
                                   err_msg='problem with time zeropoint')

class GetDateZeropointBadInput(unittest.TestCase):
    def test_get_date_zeropoint_6(self):
        test_data = MulensData()
        self.assertRaises(ValueError, test_data._get_date_zeropoint, "Potato")

    def test_get_date_zeropoint_7(self):
        test_data = MulensData()
        self.assertRaises(ValueError, test_data._get_date_zeropoint, "JD")

def test_get_jd_zeropoint_1():
    test_data = MulensData()
    assert test_data._get_jd_zeropoint(7500.) == 2450000.

def test_get_jd_zeropoint_2():
    test_data = MulensData()
    assert test_data._get_jd_zeropoint(np.array((7500., 7501.))) == 2450000.

def test_get_jd_zeropoint_3():
    test_data = MulensData()
    assert test_data._get_jd_zeropoint(np.array((57500., 57501.))) == 2400000.

def test_get_jd_zeropoint_4():
    test_data = MulensData()
    assert test_data._get_jd_zeropoint(np.array((9999., 10000.))) == 2450000.

def test_get_jd_zeropoint_5():
    test_data = MulensData()
    assert test_data._get_jd_zeropoint(np.array((2457500., 2457501.))) == 0.


class GetJDZeropointBadInput(unittest.TestCase):
    def test_get_jd_zeropoint_6(self):
        test_data = MulensData()
        self.assertRaises(ValueError, test_data._get_jd_zeropoint, 
                          np.array((np.nan, np.nan)))

    def test_get_jd_zeropoint_7(self):
        test_data = MulensData()
        self.assertRaises(ValueError, test_data._get_jd_zeropoint,
                          np.array((np.nan, 7500.)))

    def test_get_jd_zeropoint_8(self):
        test_data = MulensData()
        self.assertRaises(ValueError, test_data._get_jd_zeropoint, 
                          np.array((2450000.,7500.)))

