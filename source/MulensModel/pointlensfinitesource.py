from os.path import join, exists
import numpy as np
from scipy.interpolate import interp1d

import MulensModel as mm


class PointLensFiniteSource(object):
    """
    Data and methods used for interpolation for finite-source point-lens
    magnification calculations.

    There are no parameters of `__init__()`. In general, this class is not
    directly used by a user.
    """
    _B0B1_file_read = False

    def __init__(self):
        self._B0B1_file = join(mm.DATA_PATH, 'interpolation_table_b0b1_v3.dat')

        if not PointLensFiniteSource._B0B1_file_read:
            self._read_B0B1_file()

    def _read_B0B1_file(self):
        """Read file with pre-computed function values"""
        if not exists(self._B0B1_file):
            raise ValueError(
                'File with FSPL data does not exist.\n' + self._B0B1_file)

        kwargs = dict(usecols=range(3), unpack=True)
        (z, B0, B0_minus_B1) = np.loadtxt(self._B0B1_file, **kwargs)

        PointLensFiniteSource._B0_interpolation = interp1d(z, B0, kind='cubic')
        PointLensFiniteSource._B0_minus_B1_interpolation = interp1d(
            z, B0_minus_B1, kind='cubic')
        PointLensFiniteSource._z_min = np.min(z)
        PointLensFiniteSource._z_max = np.max(z)
        PointLensFiniteSource._B0B1_file_read = True

    def interpolate_B0(self, x):
        """
        XXX
        """
        return PointLensFiniteSource._B0_interpolation(x)

    def interpolate_B0minusB1(self, x):
        """
        XXX
        """
        return PointLensFiniteSource._B0_minus_B1_interpolation(x)

    @property
    def z_min_interpolation(self):
        """
        XXX

        *float*
        """
        return PointLensFiniteSource._z_min

    @property
    def z_max_interpolation(self):
        """
        XXX

        *float*
        """
        return PointLensFiniteSource._z_max
