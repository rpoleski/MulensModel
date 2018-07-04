import sys
import os
import unittest
import numpy as np
from astropy import units as u

import MulensModel
from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit
from MulensModel.event import Event
from MulensModel.model import Model
from MulensModel.utils import Utils


SAMPLE_FILE_01 = os.path.join(
    MulensModel.MODULE_PATH, "data", "photometry_files", "phot_ob08092_O4.dat")
SAMPLE_FILE_02 = os.path.join(
    MulensModel.MODULE_PATH, "data", "photometry_files", "OB140939",
    "ob140939_OGLE.dat")

def test_event_get_chi2_1():
    """basic unit test on ob08092 OGLE-IV data"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)
    
    ev = Event()
    mod = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data])
    ev.model = mod
    ev.datasets = [data]

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 427.20382, decimal=4, 
                                   err_msg='problem in resulting chi2')
    
    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2_no_blend), 459.09826, decimal=4, 
                                   err_msg='problem in resulting chi2 for fixed no blending')

    ev.sum_function = 'numpy.sum'
    np.testing.assert_almost_equal(ev.get_chi2(), 427.20382, decimal=4,
                                    err_msg='problem with numpy.sum')

def test_event_get_chi2_2():
    """
    Basic unit test on ob08092 OGLE-IV data. Same as above but with
    the data input twice (to test behavior for multiple datasets);
    also test if Event.get_chi2_for_dataset() gives correct answers.
    """
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    answer = 427.20382201
    
    data = MulensData(file_name=SAMPLE_FILE_01)
    
    ev = Event()
    mod = Model({'t_0':t_0, 'u_0':u_0, 't_E':t_E})
    mod.set_datasets([data, data])
    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()
    assert isinstance(chi2, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2), 2.*answer, decimal=4, 
                                   err_msg='problem in resulting chi2')
    
    chi2_no_blend = ev.get_chi2(fit_blending=False)
    assert isinstance(chi2_no_blend, float), 'wrong type of chi2'
    np.testing.assert_almost_equal(float(chi2_no_blend), 2.*459.09826, decimal=4, 
                                   err_msg='problem in resulting chi2 for fixed no blending')

    chi2_2 = ev.get_chi2_for_dataset(0)
    np.testing.assert_almost_equal(chi2_2, answer)
    
    chi2_3 = ev.get_chi2_for_dataset(1)
    np.testing.assert_almost_equal(chi2_3, answer)

def test_event_get_chi2_3():
    """test on ob08092 OGLE-IV data - MulensData.good & MulensData.bad"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    bad = np.zeros(383, dtype='bool')
    bad[300:350] = True
    data_1 = MulensData(file_name=SAMPLE_FILE_01, bad=bad)
    data_2 = MulensData(file_name=SAMPLE_FILE_01, good=~bad)
    
    ev = Event()
    mod = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    mod.set_datasets([data_1])
    ev.model = mod
    ev.datasets = [data_1]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 343.46567, decimal=4, 
                                   err_msg='problem in resulting chi2')

    mod.set_datasets([data_2])
    ev.model = mod
    ev.datasets = [data_2]
    chi2 = ev.get_chi2()
    np.testing.assert_almost_equal(float(chi2), 343.46567, decimal=4,
                                    err_msg='problem in resulting chi2')

def test_event_get_chi2_double_source_simple():
    """
    basic test on ob08092 OGLE-IV data with added second source
    Note that currently this test hacks into internal functions of 
    MulensData and MulensModel classes!
    """
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)

    ev = Event()
    mod = Model({'t_0':t_0, 'u_0':u_0, 't_E':t_E})
    mod.set_datasets([data, data])

    ev.model = mod
    ev.datasets = [data, data]

    chi2 = ev.get_chi2()
    
    assert isinstance(chi2, float), 'wrong type of chi2'
    message = 'problem in resulting chi2 for 2 exactly the same datasets'
    np.testing.assert_almost_equal(chi2, 854.407644, decimal=4, err_msg=message)

def test_event_get_chi2_3():
    """Test: If I change the model parameters, the chi2 should change."""
    #Generate a model
    t_0 = 5380.
    u_0 = 0.5
    t_E = 18.
    model = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
    
    #Generate fake data    
    times = np.arange(5320, 5420.)
    f_source = 0.1
    f_blend = 0.5
    mod_fluxes = f_source * model.magnification(times) + f_blend
    I_mag = Utils.get_mag_from_flux(mod_fluxes)
    errors = 0.01 * np.ones(times.shape)
    data = MulensData(data_list=[times, I_mag, errors])

    #Generate event and fit
    event = Event(datasets=data, model=model)
    
    orig_chi2 = event.get_chi2()
    
    #Change the model
    event.model.parameters.t_0 = 5000.

    assert event.get_chi2() != orig_chi2

def test_event_get_chi2_4():
    """test if best chi2 is remembered correctly"""
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)
    
    ev = Event()
    params = {'t_0': t_0, 'u_0': u_0, 't_E': t_E*u.day}
    mod = Model(params)
    mod.set_datasets([data])
    ev.model = mod
    ev.datasets = [data]

    chi2_1 = ev.get_chi2() # This is the best model.
    
    ev.model.parameters.parameters['t_0'] += 1.
    ev.model.parameters.parameters['u_0'] += 0.1
    ev.model.parameters.parameters['t_E'] += 1. * u.day
    chi2_2 = ev.get_chi2()

    assert chi2_2 > chi2_1
    assert ev.best_chi2 == chi2_1
    assert ev.best_chi2_parameters == params

class TestEvent(unittest.TestCase):
    def test_event_init_1(self):
        with self.assertRaises(TypeError):
            ev = Event(model=3.14)

    def test_event_init_2(self):
        with self.assertRaises(TypeError):
            ev = Event(datasets='some_string')

def test_event_chi2_gradient():
    # fs = 11.0415734, fb = 0.0 
    parameters_1 = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87}
    params_1 = ['t_0', 'u_0', 't_E']
    gradient_1 = {'t_0': 236.206598, 'u_0': 101940.249,
                't_E': -1006.88678}
    test_1 = (parameters_1, params_1, gradient_1)

    parameters_2 = {'t_0': 2456836.22, 'u_0': 0.922, 't_E': 22.87,
                    'pi_E_N': -0.248, 'pi_E_E': 0.234} # This model also
                    # used fluxes given above.
    params_2 = ['t_0', 'u_0', 't_E', 'pi_E_N', 'pi_E_E', 'f_source', 'f_blend']
    gradient_2 = {'t_0': 568.781786, 'u_0': 65235.3513, 't_E': -491.782005,
                  'pi_E_N': -187878.357, 'pi_E_E': 129162.927,
                  'f_source': -83124.5869, 'f_blend': -78653.242}
    test_2 = (parameters_2, params_2, gradient_2)
    # We're not applying the test above, yet. See 'for' loop below.

    data = MulensData(file_name=SAMPLE_FILE_02)
    kwargs = {'datasets': [data], 'coords': '17:47:12.25 -21:22:58.7'}
    
    for test in [test_1]:#, test_2]:
        (parameters, params, gradient) = test
        event = Event(model=Model(parameters), **kwargs)
        result = event.chi2_gradient(params, fit_blending=False)

        reference = np.array([gradient[key] for key in  params])
        np.testing.assert_almost_equal(reference/result, 1., decimal=1)

