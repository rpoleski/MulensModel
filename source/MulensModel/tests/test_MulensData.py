import sys
import unittest
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

from MulensModel.mulensdata import MulensData


for path in sys.path:
    if path.find("MulensModel/source") > 0:
        MODULE_PATH = "/".join(path.split("/source")[:-1])
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"


def test_file_read():
    '''read sample file and check if values match'''
    data = MulensData(file_name=SAMPLE_FILE_01, date_fmt='jdprime')

    np.testing.assert_almost_equal(data.time[0], 5264.84100, 
                                   err_msg="time of first line doesn't match")
    
    assert data.mag[-1] == 13.913, "magnitude of the last line doesn't match"

def long_test_HJD_JD_conversion():
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
    vec = np.array([1., 1.])
    dl = [vec*2455000., vec, vec*0.1]
    test_data = MulensData(data_list=dl, date_fmt="jd")
    assert test_data.time_zeropoint == 0.

def test_get_date_zeropoint_2():
    vec = np.array([1., 1.])
    dl = [vec*2455000., vec, vec*0.1]
    test_data = MulensData(data_list=dl, date_fmt="hjd")
    assert test_data.time_zeropoint == 0.

def test_get_date_zeropoint_3():
    vec = np.array([1., 1.])
    dl = [vec, vec, vec*0.1]
    test_data = MulensData(data_list=dl, date_fmt="jdprime")
    assert test_data.time_zeropoint == 2450000.

def test_get_date_zeropoint_4():
    vec = np.array([1., 1.])
    dl = [vec, vec, vec*0.1]
    test_data = MulensData(data_list=dl, date_fmt="hjdprime")
    assert test_data.time_zeropoint == 2450000.

def test_get_date_zeropoint_5():
    vec = np.array([1., 1.])
    dl = [vec*65000., vec, vec*0.1]
    test_data = MulensData(data_list=dl, date_fmt="mjd")
    assert test_data.time_zeropoint == 2400000.5
    
def test_data_list_1():
    t = np.array([7500., 7501.])
    m = np.array([21.0, 21.1])
    e = np.array([0.001, 1.000])
    data = MulensData(data_list=[t, m, e], date_fmt="jdprime")
    np.testing.assert_almost_equal(data.time, t, err_msg
                                   ='problem with time vector in MulensData')
    np.testing.assert_almost_equal(data.time_zeropoint, 2450000., 
                                   err_msg='problem with time zeropoint')

class test(unittest.TestCase):
    def test_wrong_length(self):
        with self.assertRaises(ValueError):
            t = np.array([7500., 7501.])
            m = np.array([21.0, 21.1])
            e_long = np.array([0.001, 1.000, 0.1])
            data = MulensData(data_list=[t, m, e_long], date_fmt="jdprime")

#class GetDateZeropointBadInput(unittest.TestCase):
    def test_get_date_zeropoint_6(self):
        with self.assertRaises(ValueError):
            vec = np.array([1., 1.])
            dl = [vec, vec, vec*0.1]
            test_data = MulensData(data_list=dl, date_fmt="Potato")

    def test_get_date_zeropoint_7(self):
        with self.assertRaises(ValueError):
            vec = np.array([1., 1.])
            dl = [vec, vec, vec*0.1]
            test_data = MulensData(data_list=dl, date_fmt="J_D")

