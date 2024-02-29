import os
import unittest
import numpy as np
from numpy.testing import assert_almost_equal as almost

import MulensModel as mm


dir_phot = os.path.join(mm.DATA_PATH, 'photometry_files')
SAMPLE_FILE_01 = os.path.join(dir_phot, "OB08092", "phot_ob08092_O4.dat")
SAMPLE_FILE_02 = os.path.join(dir_phot, 'OB140939', 'ob140939_Spitzer.dat')
SAMPLE_FILE_02_EPH = os.path.join(dir_phot, 'ephemeris_files',
                                  'Spitzer_ephemeris_01.dat')


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

    def test_errorbars_scaling_equation_defined(self):
        """
        make sure scaling of errorbars is done before its parameters
        are accesses
        """
        with self.assertRaises(RuntimeError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            _ = data.errorbars_scaling_equation

    def test_errorbars_scale_factors_defined(self):
        """
        make sure scaling of errorbars is done before its parameters
        are accesses
        """
        with self.assertRaises(RuntimeError):
            data = mm.MulensData(file_name=SAMPLE_FILE_01)
            _ = data.errorbars_scale_factors


def test_copy():
    """
    Test copying method
    """
    n_epochs = len(np.loadtxt(SAMPLE_FILE_01))
    random_bool = np.random.choice([False, True], n_epochs, p=[0.1, 0.9])

    data_1 = mm.MulensData(file_name=SAMPLE_FILE_01, good=random_bool)
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


def test_scale_errorbars():
    """
    Check scaling of uncertainties
    """
    factor = 1.4706
    minimum = 0.0075

    equation_1 = "sigma_mag_new = {:} * sigma_mag".format(factor)
    equation_2 = "sigma_mag_new = sqrt(sigma_mag^2 + {:}^2)".format(minimum)
    equation_3 = "sigma_mag_new = sqrt(({:} * sigma_mag)^2 + {:}^2)".format(
        factor, minimum)

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor=factor)
    almost(data.err_mag, 0.01)
    almost(data.errorbars_scale_factors, [factor, 0.])
    assert data.errorbars_scaling_equation == equation_1

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(minimum=minimum)
    almost(data.err_mag, 0.01012373)
    almost(data.errorbars_scale_factors, [1., minimum])
    assert data.errorbars_scaling_equation == equation_2

    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor, minimum)
    almost(data.err_mag, 0.0125)
    almost(data.errorbars_scale_factors, [factor, minimum])
    assert data.errorbars_scaling_equation == equation_3


def test_repr_1():
    """
    Check if one can print dataset nicely - n_bad>0
    """
    random_bad = 33 * [True] + 350 * [False]
    np.random.shuffle(random_bad)
    data = mm.MulensData(file_name=SAMPLE_FILE_01, bad=random_bad)
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5}".format(
        "phot_ob08092_O4.dat:", 383, 33)
    assert str(data) == expected


def test_repr_2():
    """
    Check if one can print dataset nicely - errorbar scaling factor
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor=1.234)
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5},".format(
        "phot_ob08092_O4.dat:", 383, 0)
    expected += " Errorbar scaling: factor = 1.234"
    assert str(data) == expected


def test_repr_3():
    """
    Check if one can print dataset nicely - bandpass and label
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_01, bandpass='I',
                         plot_properties={'label': 'OGLE'})
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5}, band = I".format(
        "OGLE:", 383, 0)
    assert str(data) == expected


def test_repr_4():
    """
    Check if one can print dataset nicely - errorbar scaling minimum
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(minimum=0.001)
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5},".format(
        "phot_ob08092_O4.dat:", 383, 0)
    expected += " Errorbar scaling: minimum = 0.001"
    assert str(data) == expected


def test_repr_5():
    """
    Check if one can print dataset nicely - ephemerides file
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_02,
                         ephemerides_file=SAMPLE_FILE_02_EPH)
    expected_begin = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5}".format(
        "ob140939_Spitzer.dat:", 31, 0)
    expected_end = "photometry_files/ephemeris_files/Spitzer_ephemeris_01.dat"
    assert str(data)[:len(expected_begin)] == expected_begin
    assert str(data)[-len(expected_end):] == expected_end


def test_repr_6():
    """
    Check if one can print dataset nicely - color
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_01,
                         plot_properties={'color': 'red'})
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5}, color = red".format(
        "phot_ob08092_O4.dat:", 383, 0)
    assert str(data) == expected


def test_repr_7():
    """
    Check if one can print dataset nicely - 2-parameter errorbar scaling
    """
    data = mm.MulensData(file_name=SAMPLE_FILE_01)
    data.scale_errorbars(factor=2.34, minimum=0.012)
    expected = "{0:25} n_epochs ={1:>5}, n_bad ={2:>5},".format(
        "phot_ob08092_O4.dat:", 383, 0)
    expected += " Errorbar scaling: factor = 2.34 minimum = 0.012"
    assert str(data) == expected
