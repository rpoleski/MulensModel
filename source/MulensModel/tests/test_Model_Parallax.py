import os
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

import MulensModel
from MulensModel.model import Model
from MulensModel.mulensdata import MulensData


DATA_PATH = os.path.join(MulensModel.MODULE_PATH, 'data')

SAMPLE_FILE_02 = os.path.join(DATA_PATH, 'phot_ob151100_OGLE_v1.dat') #HJD'
SAMPLE_FILE_02_REF = os.path.join(DATA_PATH, 'ob151100_OGLE_ref_v1.dat') #HJD'
SAMPLE_FILE_03 = os.path.join(DATA_PATH, 'phot_ob151100_Spitzer_2_v2.dat') #HJD'
SAMPLE_FILE_03_EPH = os.path.join(DATA_PATH, 'Spitzer_ephemeris_01.dat') #UTC
SAMPLE_FILE_03_REF = os.path.join(DATA_PATH, 'ob151100_Spitzer_ref_v1.dat') #HJD'

SAMPLE_ANNUAL_PARALLAX_FILE_01 = os.path.join(DATA_PATH, 'parallax_test_1.dat') #HJD'
SAMPLE_ANNUAL_PARALLAX_FILE_02 = os.path.join(DATA_PATH, 'parallax_test_2.dat') #HJD'
SAMPLE_ANNUAL_PARALLAX_FILE_03 = os.path.join(DATA_PATH, 'parallax_test_3.dat') #HJD'
SAMPLE_ANNUAL_PARALLAX_FILE_04 = os.path.join(DATA_PATH, 'parallax_test_4.dat') #HJD'
SAMPLE_ANNUAL_PARALLAX_FILE_05 = os.path.join(DATA_PATH, 'parallax_test_5.dat') #HJD'


def test_annual_parallax_calculation():
    """
    This is a high-level unit test for parallax. The "true" values were calculated from the sfit routine assuming fs=1.0, fb=0.0.
    """
    t_0 = 2457479.5 #April 1 2016, a time when parallax is large
    times = np.array([t_0-1., t_0, t_0+1., t_0+1.])
    true_no_par = [np.array([7.12399067,10.0374609, 7.12399067, 7.12399067])]
    true_with_par = [np.array([7.12376832, 10.0386009, 7.13323363, 7.13323363])]

    model_with_par = Model(t_0=t_0, u_0=0.1, t_E=10., pi_E=(0.3, 0.5),
                  coords='17:57:05 -30:22:59')
    model_with_par.parallax(satellite=False, earth_orbital=True,
                            topocentric=False)
    ones = np.ones(len(times))                    
    data = MulensData(data_list=[times, ones, ones])
    model_with_par.set_datasets([data])
    
    model_with_par.t_0_par = 2457479.
    
    model_no_par = Model(t_0=t_0, u_0=0.1, t_E=10., pi_E=(0.3, 0.5),
                  coords='17:57:05 -30:22:59')
    model_no_par.set_datasets([data])
    model_no_par.parallax(
        satellite=False, earth_orbital=False, topocentric=False)
    
    np.testing.assert_almost_equal(
        model_no_par.data_magnification, true_no_par)
    np.testing.assert_almost_equal(
        model_with_par.data_magnification, true_with_par, decimal=4)

def do_annual_parallax_test(filename):
    """testing funcations called by a few unit tests"""
    file = open(filename)
    lines = file.readlines()
    file.close()
    ulens_params = lines[3].split()
    event_params = lines[4].split()
    data = np.loadtxt(filename, dtype=None)
    model = Model(
        t_0=float(ulens_params[1])+2450000., 
        u_0=float(ulens_params[3]), 
        t_E=float(ulens_params[4]), 
        pi_E_N=float(ulens_params[5]), 
        pi_E_E=float(ulens_params[6]), 
        coords=SkyCoord(
            event_params[1]+' '+event_params[2], unit=(u.deg, u.deg)))
    model.t_0_par = float(ulens_params[2])+2450000.
    
    time = data[:,0]
    dataset = MulensData([time, 20.+time*0., 0.1+time*0.,], add_2450000=True)
    model.set_datasets([dataset])
    model.parallax(satellite=False, earth_orbital=True, topocentric=False)
    return np.testing.assert_almost_equal(
        model.data_magnification[0] / data[:,1], 1.0, decimal=4)

def test_annual_parallax_calculation_2():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_01)

def test_annual_parallax_calculation_3():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_02)

def test_annual_parallax_calculation_4():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_03)

def test_annual_parallax_calculation_5():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_04)

def test_annual_parallax_calculation_6():
    do_annual_parallax_test(SAMPLE_ANNUAL_PARALLAX_FILE_05)

#This unit test needs to be reworked so all data sets and ephemerides files are on the same 245XXXX time system.
def test_satellite_and_annual_parallax_calculation():
    """test parallax calculation with Spitzer data"""
    model_with_par = Model(t_0=2457181.93930, u_0=0.08858, t_E=20.23090, pi_E_N=-0.05413, pi_E_E=-0.16434, coords="18:17:54.74 -22:59:33.4")
    model_with_par.parallax(satellite=True, earth_orbital=True, topocentric=False)
    model_with_par.t_0_par = 2457181.9

    data_OGLE = MulensData(file_name=SAMPLE_FILE_02,add_2450000=True)
    data_Spitzer = MulensData(
        file_name=SAMPLE_FILE_03, ephemerides_file=SAMPLE_FILE_03_EPH,
        add_2450000=True)
    model_with_par.set_datasets([data_OGLE, data_Spitzer])

    ref_OGLE = np.loadtxt(SAMPLE_FILE_02_REF, unpack=True, usecols=[5])
    ref_Spitzer = np.loadtxt(SAMPLE_FILE_03_REF, unpack=True, usecols=[5])

    np.testing.assert_almost_equal(model_with_par.data_magnification[0], ref_OGLE, decimal=2)
    np.testing.assert_almost_equal(model_with_par.data_magnification[1]/ref_Spitzer, np.array([1]*len(ref_Spitzer)), decimal=3)

def test_satellite_parallax_magnification():
    """
    On a given date, the magnification should be different from the
    ground and from Spitzer. Use OB140939 as a test case and t0 as the
    test time.
    """
    t_0 = 2456836.22
    u_0 = 0.922
    t_E = 22.87
    pi_E_N = -0.248
    pi_E_E = 0.234

    ground_model = Model(t_0=t_0, u_0=u_0, t_E=t_E, pi_E=[pi_E_N, pi_E_E],
                         coords='17:47:12.25 -21:22:58.2')
    space_model =  Model(t_0=t_0, u_0=u_0, t_E=t_E, pi_E=[pi_E_N, pi_E_E], 
                         ra='17:47:12.25', dec='-21:22:58.2', 
                         ephemerides_file=SAMPLE_FILE_03_EPH)

    delta = (ground_model.magnification(t_0)
             - space_model.magnification(t_0))
    assert np.abs(delta) > 0.01
