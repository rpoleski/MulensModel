import numpy as np
from scipy import integrate
from scipy.special import ellipe
# This is incomplete elliptic integral of the second kind.

class PointLens(object):
    """
    Equations for calculating finite source effects for a point lens. 

    Takes `parameters` as input (all parameters for the model), but
    only uses rho, so maybe should modify.

    To Do:
       Add get_point_lens_magnification() method from
       MagnificationCurve. Requires also adding `methods` to the
       __init__ function.

       Might want to pass `gamma` as part of the __init__.
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
        out = 4. * z / np.pi
        function = lambda x: (1.-value**2*np.sin(x)**2)**.5

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

        function = (lambda r, theta: r * np.sqrt(1.-r**2) /
                    self.parameters.rho /
                    np.sqrt(r**2+zz**2-2.*r*zz*np.cos(theta)))
        lim_0 = lambda x: 0
        lim_1 = lambda x: 1
        W_1 = 0. * z
        for (i, zz) in enumerate(z):
            W_1[i] = integrate.dblquad(function, 0., 2.*np.pi, lim_0, lim_1)[0]

        W_1 /= np.pi
        return B_0 - 1.5 * z * self.parameters.rho * W_1

    def _get_point_lens_finite_source_magnification(
                self, u, pspl_magnification):
        """
        calculate magnification for point lens and finite source.
        The approximation was proposed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        and later the integral calculation was simplified by:

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        """
        z = u / self.parameters.rho

        magnification = pspl_magnification * self._B_0_function(z)
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
        return magnification

    def _get_point_lens_limb_darkening_magnification(
                self, u, pspl_magnification, gamma):
        """
        calculate magnification for point lens and finite source with
        limb darkening. The approximation was proposed by:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        http://adsabs.harvard.edu/abs/1994ApJ...421L..71G

        and later the integral calculation was simplified by:

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        http://adsabs.harvard.edu/abs/2004ApJ...603..139Y

        """
        z = u / self.parameters.rho
        B_0 = self._B_0_function(z)
        B_1 = self._B_1_function(z, B_0=B_0)
        magnification = pspl_magnification * (B_0 - gamma * B_1)
        return magnification

