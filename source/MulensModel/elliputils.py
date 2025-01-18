import os.path
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RegularGridInterpolator as RGI

import MulensModel as mm


class EllipUtils(object):
    """
    Data and methods used for interpolation for finite-source point-lens
    magnification calculations.

    There are no parameters of `__init__()`. In general, this class is not
    directly used by a user.
    """
    _elliptic_files_read = False

    def __init__(self):
        self.file_1_2 = os.path.join(
            mm.DATA_PATH, 'interpolate_elliptic_integral_1_2.dat')
        self.file_3 = os.path.join(
            mm.DATA_PATH, 'interpolate_elliptic_integral_3.dat')

        if not EllipUtils._elliptic_files_read:
            self._read_elliptic_files()

    def _read_elliptic_files(self):
        """
        Read 2 files with values of elliptic integrals of the 1st, 2nd,
        and 3rd kind.
        """
        (x, y1, y2) = np.loadtxt(self.file_1_2, unpack=True)
        EllipUtils._interpolate_1 = interp1d(np.log10(x), y1, kind='cubic')
        EllipUtils._interpolate_2 = interp1d(np.log10(x), y2, kind='cubic')
        EllipUtils._interpolate_1_2_x_min = np.min(np.log10(x))
        EllipUtils._interpolate_1_2_x_max = np.max(np.log10(x))

        with open(self.file_3) as file_in:
            for line in file_in.readlines():
                if line[:3] == "# X":
                    xx = np.array([float(t) for t in line.split()[2:]])
                if line[:3] == "# Y":
                    yy = np.array([float(t) for t in line.split()[2:]])

        values = np.loadtxt(self.file_3)
        try:
            EllipUtils._interpolate_3 = RGI((xx, yy), values, method='cubic', bounds_error=False)
        except ValueError:
            EllipUtils._interpolate_3 = interp2d(xx, yy, values, kind='cubic')

        EllipUtils._interpolate_3_min_x = np.min(xx)
        EllipUtils._interpolate_3_max_x = np.max(xx)
        EllipUtils._interpolate_3_min_y = np.min(yy)
        EllipUtils._interpolate_3_max_y = np.max(yy)
