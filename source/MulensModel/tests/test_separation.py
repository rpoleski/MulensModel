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
    times = get_times(parameters)

    model = mm.Model(parameters=parameters)
  
    separation = model.parameters.get_s(times)

    separation_VBB = [0.57943055, 0.36814820, 1.20000000, 1.66174048, 1.61016371]
    np.testing.assert_almost_equal(separation, separation_VBB)

 
