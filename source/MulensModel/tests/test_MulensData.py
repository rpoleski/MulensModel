import os
import unittest
import numpy as np
import warnings

import MulensModel as mm


SAMPLE_FILE_01 = os.path.join(
    mm.DATA_PATH, "photometry_files", "OB08092", "phot_ob08092_O4.dat")


def test_file_read():
    """read sample file and check if values match"""
    data = mm.MulensData(file_name=SAMPLE_FILE_01)

    np.testing.assert_almost_equal(data.time[0], 5264.84100,
                                   err_msg="time of first line doesn't match")

    assert data.mag[-1] == 13.913, "magnitude of the last line doesn't match"


def test_data_list():
    """
    check if initialization by list works fine
    """
    t = np.array([7500., 7501.])
    m = np.array([21.0, 21.1])
    e = np.array([0.001, 1.000])
    data = mm.MulensData(data_list=[t, m, e])
    np.testing.assert_almost_equal(
        data.time, t, err_msg='problem with time vector in MulensData')


class test(unittest.TestCase):
    def test_wrong_length(self):
        with self.assertRaises(ValueError):
            t = np.array([7500., 7501.])
            m = np.array([21.0, 21.1])
            e_long = np.array([0.001, 1.000, 0.1])
            _ = mm.MulensData(data_list=[t, m, e_long])

    def test_wrong_type(self):
        with self.assertRaises(TypeError):
            t = np.array([2457500., 2457501.], dtype=np.float32)
            m = np.array([21.0, 21.1])
            e = np.array([0.001, 1.000])
            _ = mm.MulensData(data_list=[t, m, e])


def test_copy():
    """
    Test copying method
    """
    n_epochs = len(np.loadtxt(SAMPLE_FILE_01))
    random_bool = np.random.choice([False, True], n_epochs, p=[0.1, 0.9])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        data_1 = mm.MulensData(file_name=SAMPLE_FILE_01, ra="18:00:00",
                               dec="-30:00:00", good=random_bool)
        data_2 = data_1.copy()

    data = [data_1.time, 100.+0.*data_1.time, 1.+0.*data_1.time]
    data_3 = mm.MulensData(data, phot_fmt='flux', bad=random_bool)
    data_4 = data_3.copy()

    assert isinstance(data_2, mm.MulensData)
    assert isinstance(data_4, mm.MulensData)

    attributes = ['time', 'mag', 'err_mag', 'flux', 'err_flux',
                  'bad', 'good', 'plot_properties']
    for attribute in attributes:
        value_1 = getattr(data_1, attribute)
        value_2 = getattr(data_2, attribute)
        assert value_1 is not value_2
        assert np.all(value_1 == value_2)
        value_1 = getattr(data_3, attribute)
        value_2 = getattr(data_4, attribute)
        assert value_1 is not value_2
        assert np.all(value_1 == value_2)

    assert data_1.coords == data_2.coords
    assert data_1.coords is not data_2.coords
    assert data_3.coords is None
    assert data_4.coords is None
