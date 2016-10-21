#! /usr/bin/env python


import sys
import unittest
import numpy as np
from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit
from MulensModel.event import Event
from MulensModel.model import Model



for path in sys.path:
    if path.find("MulensModel") > 0:
        MODULE_PATH = "/".join(path.split("/")[:-1])
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"


def test_event_get_chi2():
    '''basic unit test on ob08092 OGLE-IV data'''
    t_0 = 5379.57091
    u_0 = 0.52298
    t_E = 17.94002
    
    data = MulensData(file_name=SAMPLE_FILE_01)
    diff = (data.time - t_0) / t_E
    u2 = u_0 * u_0 + diff * diff
    u = u2**.5
    mag_model = (u2 + 2.) / (u * (u2 + 4.)**.5)
    
    ev = Event()
    mod = Model(t_0=t_0)
    mod.magnification = [mag_model]
    ev.model = mod
    ev.datasets = [data]
    chi2 = ev.get_chi2()
    #print(chi2)
    assert type(chi2) is float or type(chi2) is np.float64, 'wrong type of chi2'
    assert chi2 > 427. and chi2 < 429., 'problem in resulting chi2'
   
class TestEvent(unittest.TestCase):
    def test_event_init_1(self):
        self.assertRaises(ValueError, ev=Event(model=3.14))

    def test_event_init_2(self):
        self.assertRaises(ValueError, ev=Event(datasets='this_is_just_a_string'))

#if __name__ == "__main__":
#    test_event_get_chi2()
    
