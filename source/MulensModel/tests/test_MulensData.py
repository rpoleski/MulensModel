
import numpy as np
from MulensModel.mulensdata import MulensData


SAMPLE_FILE_01 = "../../../data/phot_ob08092_O4.dat"

def test_file_read():
    '''read sample file and check if values match'''
    data = MulensData(file_name=SAMPLE_FILE_01)

    np.testing.assert_almost_equal(data.time[0], 5264.84100, err_msg="time of first line doesn't match")
    
    assert data.mag[-1] == 13.913, "magnitude of the last line doesn't match"


