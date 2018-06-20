import numpy as np
import math
import os
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.special import ellipe
# This is an incomplete elliptic integral of the second kind.

from MulensModel.trajectory import Trajectory
import MulensModel


def get_pspl_magnification(trajectory):
    """
    This is Paczynski equation, i.e., point-source--point-lens (PSPL)
    magnification.

    Arguments :
        trajectory: *float*, *np.ndarray*, or
        *:py:class:`~MulensModel.trajectory.Trajectory`* object.
            The source-lens relative position. If _not_ a
            :py:class:`~MulensModel.trajectory.Trajectory` object,
            then trajectory is assumed to be value(s) of :math:`u`.

    Returns :
        pspl_magnification: *float* or *np.ndarray*
            The point-source--point-lens magnification for each point
            specified by `trajectory`.

    """
    if isinstance(trajectory, Trajectory):
        u2 = (trajectory.x**2 + trajectory.y**2)
    else:
        u2 = trajectory**2

    if isinstance(trajectory, float):
        pspl_magnification = (u2 + 2.) / math.sqrt(u2 * (u2 + 4.))
    else:
        pspl_magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))

    return pspl_magnification


class PointLens(object):
    """
    Equations for calculating finite source effects for a point lens.

    Keywords :
        parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
            Parameters of the model. Currently, only 
            :py:attr:`~MulensModel.modelparameters.ModelParameters.rho`
            attribute is used. 

    """

    _B0B1_file_read = False
    
    def __init__(self, parameters=None):
        self.parameters = parameters

    def _read_B0B1_file(self):
        """Read file with pre-computed function values"""
        file_ = os.path.join(MulensModel.MODULE_PATH, 'data', 
            'interpolation_table_b0b1_v1.dat')
        if not os.path.exists(file_):
            raise ValueError('File with FSPL data does not exist.\n' + file_)
        (z, B0, B0_minus_B1) = np.loadtxt(file_, unpack=True)
        PointLens._B0B1_file_read = True
        PointLens._B0_interpolation = interp1d(z, B0, kind='cubic')
        PointLens._B0_minus_B1_interpolation = interp1d(z, B0_minus_B1, 
                kind='cubic')
        PointLens._z_min = np.min(z)
        PointLens._z_max = np.max(z)
            
    def _B_0_function(self, z):
        """
        calculate B_0(z) function defined in:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        """

        out = 4. * z / np.pi
        function = lambda x: (1.-value**2*math.sin(x)**2)**.5

        for (i, value) in enumerate(z):
            if value < 1.:
                out[i] *= ellipe(value*value)
            else:
                out[i] *= integrate.quad(function, 0., np.arcsin(1./value))[0]
        return out

    def _B_1_function(self, z, B_0=None):
        """
        calculate B_1(z) function defined in:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        """
        if B_0 is None:
            B_0 = self._B_0_function(z)

        def function(r, theta):
            r_2 = r * r
            val = (1. - r_2) / (
                    r_2 + function.arg_2 + r*function.arg_3*math.cos(theta))
            return r * math.sqrt(val)

        lim_0 = lambda x: 0
        lim_1 = lambda x: 1
        rho_W_1 = 0. * z # This equals rho * W_1().
        for (i, zz) in enumerate(z):
            function.arg_1 = zz
            function.arg_2 = zz * zz
            function.arg_3 = -2. * zz
            rho_W_1[i] = integrate.dblquad(
                function, 0., 2.*np.pi, lim_0, lim_1)[0]

        rho_W_1 /= np.pi
        return B_0 - 1.5 * z * rho_W_1 

    def get_point_lens_finite_source_magnification(
                self, u, pspl_magnification, direct=False):
        """
        Calculate magnification for point lens and finite source (for
        a *uniform* source).  The approximation was proposed by:

        `Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        <http://adsabs.harvard.edu/abs/1994ApJ...421L..71G>`_

        and later the integral calculation was simplified by:

        `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        <http://adsabs.harvard.edu/abs/2004ApJ...603..139Y>`_

        Parameters :
            u: *float*, *np.array*
                The instantaneous source-lens separation. 
                Multiple values can be provided.

            pspl_magnification: *float*, *np.array*
                The point source, point lens magnification at each value of u.

            direct: *boolean*
                Use direct calculation (slow) instead of interpolation.

        Returns :
            magnification: *float*, *np.array*
                The finite source source magnification.
                Type is the same as of u parameter.

        """
        z = u / self.parameters.rho
        try:
            interator = iter(z)
        except TypeError:
            z = np.array([z])

        if not PointLens._B0B1_file_read:
            self._read_B0B1_file()

        if direct:
            mask = np.zeros_like(z, dtype=bool)
        else:
            mask = (z > PointLens._z_min) & (z < PointLens._z_max)
        
        B0 = 0. * z
        if np.any(mask): # Here we use interpolation.
            B0[mask] = PointLens._B0_interpolation(z[mask])

        mask = np.logical_not(mask)
        if np.any(mask): # Here we use direct calculation.
            B0[mask] = self._B_0_function(z[mask])

        magnification = pspl_magnification * B0
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
        return magnification

    def get_point_lens_limb_darkening_magnification(
                self, u, pspl_magnification, gamma, direct=False):
        """
        calculate magnification for point lens and finite source *with
        limb darkening*. The approximation was proposed by:

        `Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        <http://adsabs.harvard.edu/abs/1994ApJ...421L..71G>`_

        and later the integral calculation was simplified by:

        `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        <http://adsabs.harvard.edu/abs/2004ApJ...603..139Y>`_

        Parameters :
            u: *float*, *np.array*
                The instantaneous source-lens separation. Multiple values 
                can be provided.

            pspl_magnification: *float*, *np.array*
                The point source, point lens magnification at each value of u.

            gamma: *float*
                The limb-darkening coefficient. See also
                :py:class:`MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`
                
            direct: *boolean*
                Use direct calculation (very slow) instead of interpolation.

        Returns :
            magnification: *float*, *np.array*
                The finite source source magnification including
                limb-darkening. Type is the same as of u parameter.

        """
        z = u / self.parameters.rho
        try:
            interator = iter(z)
        except TypeError:
            z = np.array([z])

        if not PointLens._B0B1_file_read:
            self._read_B0B1_file()

        if direct:
            mask = np.zeros_like(z, dtype=bool)
        else:
            mask = (z > PointLens._z_min) & (z < PointLens._z_max)

        magnification = 0. * z + pspl_magnification    
        if np.any(mask): # Here we use interpolation.
            B_0 = PointLens._B0_interpolation(z[mask])
            B_0_minus_B_1 = PointLens._B0_minus_B1_interpolation(z[mask])
            magnification[mask] *= (B_0*(1.-gamma) + B_0_minus_B_1*gamma)
        
        mask = np.logical_not(mask)
        if np.any(mask): # Here we use direct calculation.
            B_0 = self._B_0_function(z[mask])
            B_1 = self._B_1_function(z[mask], B_0=B_0)
            magnification[mask] *= (B_0 - gamma * B_1)

        return magnification
