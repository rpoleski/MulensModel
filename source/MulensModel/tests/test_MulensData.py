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

    def test_scale_errorbars_twice(self):
        """make sure errorbars cannot be scaled twice"""
        with self.assertRaises(RuntimeError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            data.scale_errorbars(2.0)
            data.scale_errorbars(3.0)

    def test_scale_errorbars_negative(self):
        """make sure magnitude errobar multiplication factor is not negative"""
        with self.assertRaises(ValueError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            data.scale_errorbars(-1)

    def test_scale_errorbars_double_none(self):
        """make sure errorbar scaling gets some input"""
        with self.assertRaises(ValueError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            data.scale_errorbars()

    def test_mag_errorbars_scaling_defined(self):
        """
        make sure scaling of errorbars is done before its parameters
        are accesses
        """
        with self.assertRaises(RuntimeError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            _ = data.mag_errorbars_scaling


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


def test_scale_errorbars():
    """
    Check scaling of uncertainties
    """
    factor = 1.4706
    minimum = 0.0075

    equation_1 = "err_mag_new = factor * err_mag\n"
    equation_2 = "err_mag_new = sqrt(err_mag^2 + minimum^2)\n"
    equation_3 = "err_mag_new = sqrt((factor * err_mag)^2 + minimum^2)\n"
    factor_string = "factor: {:}".format(factor)
    minimum_string = "minimum: {:}".format(minimum)

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor=factor)
    np.testing.assert_almost_equal(data.err_mag, 0.01)
    expected_string = equation_1 + factor_string
    assert data.mag_errorbars_scaling == expected_string

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(minimum=minimum)
    np.testing.assert_almost_equal(data.err_mag, 0.01012373)
    expected_string = equation_2 + minimum_string
    assert data.mag_errorbars_scaling == expected_string

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor, minimum)
    np.testing.assert_almost_equal(data.err_mag, 0.0125)
    expected_string = equation_3 + factor_string + "\n" + minimum_string
    assert data.mag_errorbars_scaling == expected_string
