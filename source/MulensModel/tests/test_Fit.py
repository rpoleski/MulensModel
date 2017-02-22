import sys
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

import MulensModel
from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit


MODULE_PATH = "/".join(MulensModel.__file__.split("/source")[:-1])        
        
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"

def test_fit_get_input_format():
    '''read sample file and get brightness in its original format'''
    dataset = MulensData(file_name=SAMPLE_FILE_01)
    fit = Fit(data=[dataset], magnification=
                                [np.array([1.]*len(dataset.time))])                                
    input_fmt = fit.get_input_format(data=dataset)
    np.testing.assert_almost_equal(input_fmt, MulensModel.utils.MAG_ZEROPOINT, 
                     err_msg='Fit.get_input_format() returns something wrong')
   
