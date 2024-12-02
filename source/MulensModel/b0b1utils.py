from os.path import join, exists
import numpy as np
from scipy.interpolate import interp1d

import MulensModel as mm


class B0B1Utils(object):
    """
    Data and methods used for interpolation for finite-source point-lens
    magnification calculations.

    There are no parameters of `__init__()`. In general, this class is not
    directly used by a user.
    """
    _B0B1_file_read = False

    def __init__(self):
        self._B0B1_file = join(mm.DATA_PATH, 'interpolation_table_b0b1_v3.dat')

        if not B0B1Utils._B0B1_file_read:
            self._read_B0B1_file()

    def _read_B0B1_file(self):
        """Read file with pre-computed function values"""
        if not exists(self._B0B1_file):
            raise ValueError(
                'File with FSPL data does not exist.\n' + self._B0B1_file)

        file_info = np.loadtxt(self._B0B1_file, unpack=True)
        (z, B0, B0_minus_B1, B1, B0_prime, B1_prime) = file_info

        kwargs = {'kind': 'cubic', 'bounds_error': False, 'fill_value': 1.}
        B0B1Utils._B0_interpolation = interp1d(z, B0, **kwargs)
        kwargs['fill_value'] = 0.
        B0B1Utils._B0_minus_B1_interpolation = interp1d(
            z, B0_minus_B1, **kwargs)
        B0B1Utils._B1_interpolation = interp1d(z, B1, **kwargs)
        B0B1Utils._B0_prime_interpolation = interp1d(
            z, B0_prime, **kwargs)
        B0B1Utils._B1_prime_interpolation = interp1d(
            z, B1_prime, **kwargs)
        B0B1Utils._z_min = np.min(z)
        B0B1Utils._z_max = np.max(z)
        B0B1Utils._B0B1_file_read = True

    def interpolate_B0(self, z):
        """
        Interpolate B_0(z) function.

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            B0: *np.ndarray*
                Values of B_0.
        """
        return B0B1Utils._B0_interpolation(z)

    def interpolate_B0minusB1(self, z):
        """
        Interpolate B_0(z)-B_1(z) function.

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            B0_minus_B1: *np.ndarray*
                Values of B_0(z)-B_1(z).
        """
        return B0B1Utils._B0_minus_B1_interpolation(z)

    def interpolate_B1(self, z):
        """
        Interpolate B_1(z) function.

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            B1: *np.ndarray*
                Values of B_1(z).
        """
        return B0B1Utils._B1_interpolation(z)

    def interpolate_B0prime(self, z):
        """
        Interpolate derivative of B_0(z), i.e., d B_0(z)/d z.

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            dB0_dz: *np.ndarray*
                Values of d B_0(z)/d z.
        """
        return B0B1Utils._B0_prime_interpolation(z)

    def interpolate_B1prime(self, z):
        """
        Interpolate derivative of B_1(z), i.e., d B_1(z)/d z.

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            dB1_dz: *np.ndarray*
                Values of d B_1(z)/d z.
        """
        return B0B1Utils._B1_prime_interpolation(z)

    def get_interpolation_mask(self, z):
        """
        Get mask of z values for which interpolation of B_0(z), B_1(z) etc.
        is allowed

        Parameters :
            z: *np.ndarray*
                Values of u/rho.

        Returns :
            mask: *np.ndarray*
                Mask indicating for which z values interpolation is allowed.
        """
        mask = (z > B0B1Utils._z_min)
        mask &= (z < B0B1Utils._z_max)
        return mask

    @property
    def z_max_interpolation(self):
        """
        Maximum value of z for which interpolation is allowed.

        *float*
        """
        return B0B1Utils._z_max
