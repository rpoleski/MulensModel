import MulensModel as mm
import numpy as np
import unittest
import pytest
from test_compare_trajectory import get_parameters, get_times

def test_separation():
    """
    compares separation to values form VBBinaryLensing
    """
    parameters = get_parameters()
    times = get_times()

    model = mm.Model(parameters=parameters)
  
    separation = model.parameters.get_s(times)
    
    separation_VBB = [0.918651, 1.332512, 1.200000, 1.079128, 1.427853]
    np.testing.assert_almost_equal(separation, separation_VBB)

 
