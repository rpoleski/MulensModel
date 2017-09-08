import sys, os
import unittest
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

import MulensModel
from MulensModel.mulensdata import MulensData

        
SAMPLE_FILE_01 = os.path.join(MulensModel.MODULE_PATH, 
                                    "data", "phot_ob08092_O4.dat")

def test_file_read():
    """read sample file and check if values match"""
    data = MulensData(file_name=SAMPLE_FILE_01)

    np.testing.assert_almost_equal(data.time[0], 5264.84100, 
                                   err_msg="time of first line doesn't match")
    
    assert data.mag[-1] == 13.913, "magnitude of the last line doesn't match"
 
def test_data_list_1():
    t = np.array([7500., 7501.])
    m = np.array([21.0, 21.1])
    e = np.array([0.001, 1.000])
    data = MulensData(data_list=[t, m, e])
    np.testing.assert_almost_equal(data.time, t, err_msg
                                   ='problem with time vector in MulensData')

class test(unittest.TestCase):
    def test_wrong_length(self):
        with self.assertRaises(ValueError):
            t = np.array([7500., 7501.])
            m = np.array([21.0, 21.1])
            e_long = np.array([0.001, 1.000, 0.1])
            data = MulensData(data_list=[t, m, e_long])

