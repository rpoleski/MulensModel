import unittest
from MulensModel.mulenstime import MulensTime

def test_get_date_zeropoint_1():
    test_data = MulensTime()
    assert test_data._get_date_zeropoint(date_fmt="jd") == 0.

def test_get_date_zeropoint_2():
    test_data = MulensTime()
    assert test_data._get_date_zeropoint(date_fmt="hjd") == 0.

def test_get_date_zeropoint_3():
    test_data = MulensTime()
    assert test_data._get_date_zeropoint(date_fmt="jdprime") == 2450000.

def test_get_date_zeropoint_4():
    test_data = MulensTime()
    assert test_data._get_date_zeropoint(date_fmt="hjdprime") == 2450000.

def test_get_date_zeropoint_5():
    test_data = MulensTime()
    assert test_data._get_date_zeropoint(date_fmt="mjd") == 2400000.5

class GetDateZeropointBadInput(unittest.TestCase):
    def test_get_date_zeropoint_6(self):
        test_data = MulensTime()
        self.assertRaises(ValueError, test_data._get_date_zeropoint, "Potato")

    def test_get_date_zeropoint_7(self):
        test_data = MulensTime()
        self.assertRaises(ValueError, test_data._get_date_zeropoint, "JD")

