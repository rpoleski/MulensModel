import numpy as np
from math import sin, cos, sqrt
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
        :py:class:`~MulensModel.trajectory.Trajectory` object

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
        pspl_magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
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
        file_ = os.path.join(
            MulensModel.MODULE_PATH, 'data',
            'interpolation_table_b0b1_v1.dat')
        if not os.path.exists(file_):
            file_ = os.path.join(
                os.path.dirname(__file__), 'data',
                'interpolation_table_b0b1_v1.dat')
            if not os.path.exists(file_):
                raise ValueError(
                    'File with FSPL data does not exist.\n' + file_)
        (z, B0, B0_minus_B1) = np.loadtxt(file_, unpack=True)
        PointLens._B0B1_file_read = True
        PointLens._B0_interpolation = interp1d(z, B0, kind='cubic')
        PointLens._B0_minus_B1_interpolation = interp1d(
                z, B0_minus_B1, kind='cubic')
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
        function = lambda x: (1.-value**2*sin(x)**2)**.5

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
                    r_2 + function.arg_2 + r*function.arg_3*cos(theta))
            return r * sqrt(val)

        lim_0 = lambda x: 0
        lim_1 = lambda x: 1
        rho_W_1 = 0. * z  # This equals rho * W_1().
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
            _ = iter(z)
        except TypeError:
            z = np.array([z])

        if not PointLens._B0B1_file_read:
            self._read_B0B1_file()

        if direct:
            mask = np.zeros_like(z, dtype=bool)
        else:
            mask = (z > PointLens._z_min) & (z < PointLens._z_max)

        B0 = 0. * z
        if np.any(mask):  # Here we use interpolation.
            B0[mask] = PointLens._B0_interpolation(z[mask])

        mask = np.logical_not(mask)
        if np.any(mask):  # Here we use direct calculation.
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
                :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`

            direct: *boolean*
                Use direct calculation (very slow) instead of interpolation.

        Returns :
            magnification: *float*, *np.array*
                The finite source source magnification including
                limb-darkening. Type is the same as of u parameter.

        """
        z = u / self.parameters.rho
        try:
            _ = iter(z)
        except TypeError:
            z = np.array([z])

        if not PointLens._B0B1_file_read:
            self._read_B0B1_file()

        if direct:
            mask = np.zeros_like(z, dtype=bool)
        else:
            mask = (z > PointLens._z_min) & (z < PointLens._z_max)

        magnification = 0. * z + pspl_magnification
        if np.any(mask):  # Here we use interpolation.
            B_0 = PointLens._B0_interpolation(z[mask])
            B_0_minus_B_1 = PointLens._B0_minus_B1_interpolation(z[mask])
            magnification[mask] *= (B_0*(1.-gamma) + B_0_minus_B_1*gamma)

        mask = np.logical_not(mask)
        if np.any(mask):  # Here we use direct calculation.
            B_0 = self._B_0_function(z[mask])
            B_1 = self._B_1_function(z[mask], B_0=B_0)
            magnification[mask] *= (B_0 - gamma * B_1)

        return magnification

    def get_point_lens_uniform_integrated_magnification(self, u, rho):
        """
        XXX

        `Lee, C.-H. et al. 2009 ApJ 695, 200 "Finite-Source Effects in
        Microlensing: A Precise, Easy to Implement, Fast, and Numerically
        Stable Formalism"
        <http://adsabs.harvard.edu/abs/2009ApJ...695..200L>`_
        """
        n = 100

        mag = np.zeros_like(u)

        for i in range(len(u)):
            if u[i] > rho:
                mag[i] = self._temp_noLD_2a(u[i], rho, n)
            else:
                mag[i] = self._temp_noLD_2b(u[i], rho, n)
        return mag

    def _temp_noLD_2a(self, u, rho, n):
        """XXX"""
        if n % 2 != 0:
            raise ValueError('odd number expected')
        theta_max = np.arcsin(rho / u)
        rho2 = rho * rho

        def u_1(theta):
            """XXX here theta is a vector; Eq. 4 of Lee+09"""
            if u <= rho:
                return 0.
            out = np.zeros_like(theta)
            mask = (theta <= theta_max)
            if np.any(mask):
                ucos = u * np.cos(theta[mask])
                out[mask] = ucos - np.sqrt(rho2 - u * u + ucos**2)
            return out

        def u_2(theta):
            """XXX here theta is a vector; Eq. 5 of Lee+09"""
            if u <= rho:
                ucos = u * np.cos(theta)
                return ucos + np.sqrt(rho2 - u * u + ucos**2)
            else:
                out = np.zeros_like(theta)
                mask = (theta <= theta_max)
                if np.any(mask):
                    ucos = u * np.cos(theta[mask])
                    out[mask] = ucos + np.sqrt(rho2 - u * u + ucos**2)
                return out

        def f(theta):
            """
            XXX
            here theta is a vector
            equation in text between Eq. 7 and 8
            """
            u_1_ = u_1(theta)
            u_2_ = u_2(theta)
            return u_2_*np.sqrt(u_2_**2+4.) - u_1_*np.sqrt(u_1_**2+4.)

        out = (u+rho)*sqrt((u+rho)**2+4.)-(u-rho)*sqrt((u-rho)**2+4.)
        vector_1 = np.arange(1., (n/2 - 1.) + 1)
        vector_2 = np.arange(1., n/2 + 1)
        out += 2. * np.sum(f(2.*vector_1*theta_max/n))
        out += 4. * np.sum(f((2.*vector_2-1.)*theta_max/n))
        out *= theta_max / (3.*np.pi*rho2*n)
        return out

    def _temp_noLD_2b(self, u, rho, n):
        """XXX"""
        if n % 2 != 0:
            raise ValueError('odd number expected')
        rho2 = rho * rho

        def u_1(theta):
            """XXX here theta is a vector; Eq. 4 of Lee+09"""
            if u <= rho:
                return 0.
            out = np.zeros_like(theta)
            mask = (theta <= theta_max)
            if np.any(mask):
                ucos = u * np.cos(theta[mask])
                out[mask] = ucos - np.sqrt(rho2 - u * u + ucos**2)
            return out

        def u_2(theta):
            """XXX here theta is a vector; Eq. 5 of Lee+09"""
            if u <= rho:
                ucos = u * np.cos(theta)
                return ucos + np.sqrt(rho2 - u * u + ucos**2)
            else:
                out = np.zeros_like(theta)
                mask = (theta <= theta_max)
                if np.any(mask):
                    ucos = u * np.cos(theta[mask])
                    out[mask] = ucos + np.sqrt(rho2 - u * u + ucos**2)
                return out

        def f(theta):
            """
            XXX
            here theta is a vector;
            equation in text between Eq. 7 and 8
            """
            u_1_ = u_1(theta)
            u_2_ = u_2(theta)
            return u_2_*np.sqrt(u_2_**2+4.) - u_1_*np.sqrt(u_1_**2+4.)

        out = (u+rho)*sqrt((u+rho)**2+4.)-(u-rho)*sqrt((u-rho)**2+4.)
        vector_1 = np.arange(1., (n - 1.) + 1)
        vector_2 = np.arange(1., n + 1)
        out += 2. * np.sum(f(vector_1 * np.pi / n))
        out += 4. * np.sum(f((2.*vector_2-1.)*np.pi/(2.*n)))
        out /= 2. * 3. * n * rho2
        return out

####################################################
# OBSOLETE CODE BELOW
####################################################

    def _temp_noLD(self, u, rho):
        """XXX"""
        if u >= rho:
            theta_max = np.arcsin(rho / u)
        rho2 = rho * rho
        u2 = u * u

        def u_1(theta):
            """XXX"""
            if u <= rho:
                return 0.
            if theta > theta_max:
                return 0.
            return u*cos(theta) - sqrt(rho2 - (u*sin(theta))**2)
#            return u_1.out

        def u_2(theta):
            """XXX"""
            if u > rho and theta > theta_max:
                return 0.
            return u*cos(theta) + sqrt(rho2 - (u*sin(theta))**2)
#            a = u*cos(theta)
#            b = sqrt(rho2 - u2 + a * a)
#            u_1.out = a - b
#            return a + b

        def short_fun(u_):
            """XXX"""
            if u_ == 0.:
                return 0.
            return u_ * sqrt(u_ * u_ + 4.)

        def integrand(theta):
            """XXX"""
            u_1_ = u_1(theta)
            u_2_ = u_2(theta)
            return short_fun(u_2_) - short_fun(u_1_)
#            return u_2_*sqrt(u_2_**2+4.) - u_1_*sqrt(u_1_**2+4.)

        integral = integrate.quad(integrand, 0, np.pi, epsabs=0., epsrel=1e-5)
        return integral[0] / (np.pi * rho2)



    def _temp(self, u, rho, gamma):
        """
        apply general method from Lee+09 to a single datapoint
        """
        if u >= rho:
            theta_max = np.arcsin(rho / u)
        rho2 = rho * rho
        u2 = u * u

        def u_1(theta):
            """XXX"""
            if u <= rho:
                return 0.
            if theta > theta_max:
                return 0.
#            return u_1.out
            return u*cos(theta) - sqrt(rho2 - (u*sin(theta))**2)

        def u_2(theta):
            """XXX"""
            if u > rho and theta > theta_max:
                return 0.
#            a = u*cos(theta)
#            b = sqrt(rho2 - u2 + a * a)
#            u_1.out = a - b
#            return a + b
            return u*cos(theta) + sqrt(rho2 - (u*sin(theta))**2)

        def temp(u_2):
            return (u_2 + 2.) / sqrt(u_2 + 4.)

#        fac = 1.5 / rho

        def fun(u_, theta):
            """XXX"""
#            if fun.theta != theta:
#                fun.theta = theta
#                fun.temp = 2. * cos(theta)
            u_2 = u_ * u_
            frac = (u_2 - 2.*u*u_*cos(theta) + u2) / rho2
#            frac = u_2 - 2.*u*u_*cos(theta) + u2
#            if frac > rho2:
#                return 0.
#            return temp(u_2)*(1.-gamma*(1.-fac*sqrt(rho2-frac)))
#            frac = (u_2 - fun.temp*u*u_ + u2) / rho2
            if frac > 1.:
                return 0.
#            return temp(u_2)*(1.-gamma*(1.-1.5*sqrt(1.-frac)))
            temp = (u_2 + 2.) / sqrt(u_2 + 4.)
            return temp*(1.-gamma*(1.-1.5*sqrt(1.-frac)))

#        fun.theta = None

        e = 1e-2
        integral = integrate.dblquad(fun, 0., np.pi, u_1, u_2,
                                     epsabs=0., epsrel=e)
        mag = integral[0] * 2. / (np.pi * rho * rho)
        return mag

