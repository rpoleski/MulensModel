import numpy as np
import math
from scipy import integrate
from scipy.special import ellipe
# This is an incomplete elliptic integral of the second kind.

from MulensModel.trajectory import Trajectory


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

    def __init__(self, parameters=None):
        self.parameters = parameters

    def _B_0_function(self, z):
        """
        calculate B_0(z) function defined in:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        """
        try:
            interator = iter(z)
        except TypeError:
            z = np.array([z])

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
        try:
            interator = iter(z)
        except TypeError:
            z = np.array([z])

        if B_0 is None:
            B_0 = self._B_0_function(z)

        def function(r, theta):
            r_2 = r * r
            val = (1. - r_2) / (
                    r_2 + function.arg_2 + r*function.arg_3*math.cos(theta))
            return r * function.arg_4 * math.sqrt(val)

        lim_0 = lambda x: 0
        lim_1 = lambda x: 1
        W_1 = 0. * z
        for (i, zz) in enumerate(z):
            function.arg_1 = zz
            function.arg_2 = zz * zz
            function.arg_3 = -2. * zz
            function.arg_4 = 1. / self.parameters.rho
            W_1[i] = integrate.dblquad(function, 0., 2.*np.pi, lim_0, lim_1)[0]

        W_1 /= np.pi
        return B_0 - 1.5 * z * self.parameters.rho * W_1

    def get_point_lens_finite_source_magnification(
                self, u, pspl_magnification):
        """
        Calculate magnification for point lens and finite source (for
        a *uniform* source).  The approximation was proposed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        and later the integral calculation was simplified by:

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        Parameters :
            u: *float*, *np.array*
                The instantaneous source-lens separation. 
                Multiple values can be provided.

            pspl_magnification: *float*, *np.array*
                The point source, point lens magnification at each value of u.

        Returns :
            magnification: *float*, *np.array*
                The finite source source magnification.
                Type is the same as of u parameter.

        """
        z = u / self.parameters.rho

        magnification = pspl_magnification * self._B_0_function(z)
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
        return magnification

    def get_point_lens_limb_darkening_magnification(
                self, u, pspl_magnification, gamma):
        """
        calculate magnification for point lens and finite source *with
        limb darkening*. The approximation was proposed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        and later the integral calculation was simplified by:

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        Parameters :
            u: *float*, *np.array*
                The instantaneous source-lens separation. Multiple values 
                can be provided.

            pspl_magnification: *float*, *np.array*
                The point source, point lens magnification at each value of u.

            gamma: *float*
                The limb-darkening coefficient. See also
                :py:class:`MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`

        Returns :
            magnification: *float*, *np.array*
                The finite source source magnification including
                limb-darkening. Type is the same as of u parameter.

        """
        z = u / self.parameters.rho
        B_0 = self._B_0_function(z)
        B_1 = self._B_1_function(z, B_0=B_0)
        magnification = pspl_magnification * (B_0 - gamma * B_1)
        return magnification
