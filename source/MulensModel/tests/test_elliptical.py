import unittest
import pytest
import numpy as np

import MulensModel as mm


def setup_keplerian_elliptical(dict_to_add):
    """
    Setup dictionary for tests of elliptical keplerian motion
    """
    dict_static = {'t_0': 2456789.01234, 'u_0': 1., 't_E': 12.345,
                   's': 1.2345, 'q': 0.01234, 'alpha': 30., 'rho': 0.001, 'ds_dt': 0.1, 'dalpha_dt': 10. }

    return {**dict_static, **dict_to_add}


class test_Keplerian_elliptical(unittest.TestCase):
    def test_keplerian_no_s_z(self):
        """fails if s_z is not given"""
        dict_2 = setup_keplerian_elliptical({'ds_z_dt': 1.9, 'a_r': 1.22})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_2)

    def test_keplerian_no_ds_z_dt(self):
        """fails if ds_z_dt is not given"""
        dict_3 = setup_keplerian_elliptical({'s_z': 0.1, 'a_r': 1.22})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_3)

    def test_keplerian_only_z(self):
        """fails if s_z and ds_z_dt are given only"""
        dict_4 = setup_keplerian_elliptical({'s_z': 0.1, 'ds_z_dt': 1.9})
        with self.assertRaises(KeyError):
            mm.ModelParameters(dict_4)

