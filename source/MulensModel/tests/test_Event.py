
import sys
import unittest
import numpy as np

import MulensModel
from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit
from MulensModel.event import Event
from MulensModel.model import Model
from MulensModel.utils import Utils


MODULE_PATH = "/".join(MulensModel.__file__.split("/source")[:-1])
        
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"


def test_event_get_chi2():
    '''basic unit test on ob08092 OGLE-IV data'''
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)
    
    ev = Event()
    mod = Model(t_0=t_0, u_0=u_0, t_E=t_E)
    mod.set_datasets([data])
    ev.model = mod
    ev.datasets = [data]

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 428.58655, decimal=4, 
                                   err_msg='problem in resulting chi2')
    
    chi2_no_blend = ev.get_chi2(fit_blending_all=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2_no_blend), 460.72308, decimal=4, 
                                   err_msg='problem in resulting chi2 for fixed no blending')


def test_event_get_chi2_double_source_simple():
    '''basic test on ob08092 OGLE-IV data with added second source
    Note that currently this test hacks into internal functions of 
    MulensData and MulensModel classes!'''
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)

    ev = Event()
    mod = Model(t_0=t_0, u_0=u_0, t_E=t_E)
    mod.set_datasets([data, data])

    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()
    
    assert isinstance(chi2, float), 'wrong type of chi2'
    message = 'problem in resulting chi2 for 2 exactly the same datasets'
    np.testing.assert_almost_equal(chi2, 857.17310, decimal=4, err_msg=message)

"""
def test_event_get_chi2():
    #Generate a model
    t_0 = 5380.
    u_0 = 0.5
    t_E = 18.
    model = Model(t_0=t_0, u_0=u_0, t_E=t_E)
    
    #Generate fake data offset from that model
    times = np.arange(5320,5420.)

    f_source = 0.1
    f_blend = 0.5
    mod_fluxes = f_source * model.magnification(times) + f_blend

    I_mag = Utils.get_mag_from_flux(mod_fluxes)

    random_numbers = np.random.randn(times.size)
    errors = 0.01 * np.ones(times.shape)
    
    mags = I_mag + random_numbers * errors

    data = MulensData(data_list=[times, mags, errors])

    #Generate event and fit
    event = Event(datasets=data, model=model)
    
    print('This unit test does not work. Not sure if it is designed correctly.')

    expected_chi2 = np.sum(np.abs(random_numbers)**2

    #assert np.abs(( event.get_chi2() - expected_chi2 ) / expected_chi2) < 0.05

    #Change the model
    event.model.t_0 = 5000.

    data_fluxes, data_flux_errors = Utils.get_flux_and_err_from_mag(mags,errors)
    mean_flux = np.mean(data_fluxes)
    chi2_expected = 0.
    for i, flux in enumerate(data_fluxes):
        chi2_expected += ( (flux - mean_flux) / data_flux_errors[i] )**2

    #assert np.abs(( event.get_chi2() - chi2_expected ) / chi2_expected) < 0.05
"""

class TestEvent(unittest.TestCase):
    def test_event_init_1(self):
        with self.assertRaises(TypeError):
            ev = Event(model=3.14)

    def test_event_init_2(self):
        with self.assertRaises(TypeError):
            ev = Event(datasets='some_string')

