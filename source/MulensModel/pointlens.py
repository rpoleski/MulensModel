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
        https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract

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
        https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract

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
        <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_

        and later the integral calculation was simplified by:

        `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        <https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract>`_

        This approach assumes rho is small (rho < 0.1). For larger sources
        use :py:func:`get_point_lens_uniform_integrated_magnification`.

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
        <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_

        and later the integral calculation was simplified by:

        `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        <https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract>`_

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
        Calculate magnification for the point lens and uniform finite source.
        This approach works well for for small and large sources
        (e.g., rho~0.5). Uses the method presented by:

        `Lee, C.-H. et al. 2009 ApJ 695, 200 "Finite-Source Effects in
        Microlensing: A Precise, Easy to Implement, Fast, and Numerically
        Stable Formalism"
        <https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract>`_
        """
        n = 100

        mag = np.zeros_like(u)

        for i in range(len(u)):
            if u[i] > rho:
                mag[i] = self._noLD_Lee09_large_u(u[i], rho, n)
            else:
                mag[i] = self._noLD_Lee09_small_u(u[i], rho, n)
        return mag

    def _u_1_Lee09(self, theta, u, rho, theta_max=None):
        """
        Calculates Equation 4 of Lee et al. 2009.
        The u variable is float, theta is np.ndarray.
        """
        if u <= rho:
            return 0. * theta
        out = np.zeros_like(theta)
        mask = (theta <= theta_max)
        if np.any(mask):
            ucos = u * np.cos(theta[mask])
            out[mask] = ucos - np.sqrt(rho * rho - u * u + ucos**2)
        return out

    def _u_2_Lee09(self, theta, u, rho, theta_max=None):
        """
        Calculates Equation 5 of Lee et al. 2009.
        The u variable is float, theta is np.ndarray.
        """
        if u <= rho:
            ucos = u * np.cos(theta)
            return ucos + np.sqrt(rho * rho - u * u + ucos**2)
        else:
            out = np.zeros_like(theta)
            mask = (theta <= theta_max)
            if np.any(mask):
                ucos = u * np.cos(theta[mask])
                out[mask] = ucos + np.sqrt(rho * rho - u * u + ucos**2)
            return out

    def _f_Lee09(self, theta, u, rho, theta_max=None):
        """
        Calculates equation in text between Eq. 7 and 8 from
        Lee et al. 2009.
        """
        u_1_ = self._u_1_Lee09(theta, u, rho, theta_max)
        u_2_ = self._u_2_Lee09(theta, u, rho, theta_max)

        f_u_1 = u_1_ * np.sqrt(u_1_**2 + 4.)
        f_u_2 = u_2_ * np.sqrt(u_2_**2 + 4.)
        return f_u_2 - f_u_1

    def _noLD_Lee09_large_u(self, u, rho, n):
        """
        Calculates Equation 7 from Lee et al. 2009 in case u > rho.
        """
        if n % 2 != 0:
            raise ValueError('internal error - odd number expected')
        theta_max = np.arcsin(rho / u)
        out = (u+rho)*sqrt((u+rho)**2+4.)-(u-rho)*sqrt((u-rho)**2+4.)
        vector_1 = np.arange(1., (n/2 - 1.) + 1)
        vector_2 = np.arange(1., n/2 + 1)
        arg_1 = 2. * vector_1 * theta_max / n
        arg_2 = (2. * vector_2 - 1.) * theta_max / n
        out += 2. * np.sum(self._f_Lee09(arg_1, u, rho, theta_max))
        out += 4. * np.sum(self._f_Lee09(arg_2, u, rho, theta_max))
        out *= theta_max / (3. * np.pi * rho * rho * n)
        return out

    def _noLD_Lee09_small_u(self, u, rho, n):
        """
        Calculates Equation 7 from Lee et al. 2009 in case u < rho.
        """
        if n % 2 != 0:
            raise ValueError('internal error - odd number expected')
        out = (u+rho)*sqrt((u+rho)**2+4.)-(u-rho)*sqrt((u-rho)**2+4.)
        vector_1 = np.arange(1., (n - 1.) + 1)
        vector_2 = np.arange(1., n + 1)
        arg_1 = vector_1 * np.pi / n
        arg_2 = (2. * vector_2 - 1.) * np.pi / (2. * n)
        out += 2. * np.sum(self._f_Lee09(arg_1, u, rho))
        out += 4. * np.sum(self._f_Lee09(arg_2, u, rho))
        out /= 2. * 3. * n * rho * rho
        return out

    def get_point_lens_LD_integrated_magnification(self, u, rho, gamma):
        """
        Calculate magnification for the point lens and finite source with
        limb-darkening. This approach works well for for small and large
        sources (e.g., rho~0.5). Uses the method presented by:

        `Lee, C.-H. et al. 2009 ApJ 695, 200 "Finite-Source Effects in
        Microlensing: A Precise, Easy to Implement, Fast, and Numerically
        Stable Formalism"
        <https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract>`_
        """
        n_theta = 90
        n_u = 1000

        mag = np.zeros_like(u)

        for i in range(len(u)):
            mag[i] = self._LD_Lee09(u[i], rho, gamma, n_theta, n_u)

        return mag

    def _integrand_Lee09_v2(self, u_, u, theta_, rho, gamma):
        """
        Integrand in Equation 13 in Lee et al. 2009.

        u_ and theta_ are np.ndarray, other parameters are scalars.
        theta_ is in fact cos(theta_) here
        """
        values = 1. - (u_**2 - 2. * u_ * u * theta_ + u**2) / rho**2
        values[:, -1] = 0.
        if values[-1, 0] < 0.:
            values[-1, 0] = 0.
        out = 1. - gamma * (1. - 1.5 * np.sqrt(values))
        return out * (u_**2 + 2.) / np.sqrt(u_**2 + 4.)

    def _LD_Lee09(self, u, rho, gamma, n_theta, n_u):
        """
        Calculates Equation 13 from Lee et al. 2009.
        """
        accuracy = 1e-4
        theta_sub = 1.e-12
        u_1_min = 1.e-13

        if n_theta % 2 != 0:
            raise ValueError('internal error - odd number expected')
        if n_u % 2 != 0:
            raise ValueError('internal error - odd number expected')

        if u > rho:
            theta_max = np.arcsin(rho / u)
        else:
            theta_max = np.pi
        theta = np.linspace(0, theta_max-theta_sub, n_theta)
        integrand_values = np.zeros_like(theta)
        u_1 = self._u_1_Lee09(theta, u, rho, theta_max)
        u_1 += u_1_min
        u_2 = self._u_2_Lee09(theta, u, rho, theta_max)

        size = (len(theta), n_u)
        temp = np.zeros(size)
        temp2 = (np.zeros(size).T + np.cos(theta)).T
        for (i, (theta_, u_1_, u_2_)) in enumerate(zip(theta, u_1, u_2)):
            temp[i] = np.linspace(u_1_, u_2_, n_u)
        integrand = self._integrand_Lee09_v2(temp, u, temp2, rho, gamma)
        for (i, (theta_, u_1_, u_2_)) in enumerate(zip(theta, u_1, u_2)):
            integrand_values[i] = integrate.simps(integrand[i],
                                                  dx=temp[i, 1]-temp[i, 0])
        out = integrate.simps(integrand_values, dx=theta[1]-theta[0])
        out *= 2. / (np.pi * rho**2)
        return out

