import os
import warnings
import numpy as np
from math import sin, cos, sqrt, log10
from scipy import integrate
from scipy.interpolate import interp1d, interp2d
from scipy.special import ellipk, ellipe
# These are complete elliptic integrals of the first and the second kind.
from sympy.functions.special.elliptic_integrals import elliptic_pi as ellip3

import MulensModel as mm


# def get_pspl_magnification(trajectory):
#     """
#     This is Paczynski equation, i.e., point-source--point-lens (PSPL)
#     magnification.
#
#     Arguments :
#         trajectory: *float*, *np.ndarray*, or
#         :py:class:`~MulensModel.trajectory.Trajectory` object
#
#             The source-lens relative position. If _not_ a
#             :py:class:`~MulensModel.trajectory.Trajectory` object,
#             then trajectory is assumed to be value(s) of :math:`u`.
#
#     Returns :
#         pspl_magnification: *float* or *np.ndarray*
#             The point-source--point-lens magnification for each point
#             specified by `trajectory`.
#
#     """
#     if isinstance(trajectory, mm.Trajectory):
#         u2 = (trajectory.x**2 + trajectory.y**2)
#     else:
#         u2 = trajectory**2
#
#     if isinstance(trajectory, float):
#         pspl_magnification = (u2 + 2.) / sqrt(u2 * (u2 + 4.))
#     else:
#         pspl_magnification = (u2 + 2.) / np.sqrt(u2 * (u2 + 4.))
#
#     return pspl_magnification
#
#
class PointLens(object):
#     """
#     Equations for calculating finite source effects for a point lens.
#
#     Keywords :
#         parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
#             Parameters of the model. Currently, only
#             :py:attr:`~MulensModel.modelparameters.ModelParameters.rho`
#             attribute is used.
#
#     Attributes :
#         parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
#             input parameters
#     """
#
#     _elliptic_files_read = False
#
     def __init__(self, parameters=None):
         pass
#         if not isinstance(parameters, mm.ModelParameters):
#             raise TypeError(
#                 "PointLens argument has to be of ModelParameters type, not " +
#                 str(type(parameters)))
#
#         self.parameters = parameters
#         self._B0B1_data = mm.B0B1Utils()
#
#     def _read_elliptic_files(self):
#         """
#         Read 2 files with values of elliptic integrals of the 1st, 2nd,
#         and 3rd kind.
#         """
#         file_1_2 = os.path.join(
#             mm.DATA_PATH, 'interpolate_elliptic_integral_1_2.dat')
#         file_3 = os.path.join(
#             mm.DATA_PATH, 'interpolate_elliptic_integral_3.dat')
#
#         (x, y1, y2) = np.loadtxt(file_1_2, unpack=True)
#         PointLens._interpolate_1 = interp1d(np.log10(x), y1, kind='cubic')
#         PointLens._interpolate_2 = interp1d(np.log10(x), y2, kind='cubic')
#         PointLens._interpolate_1_2_x_min = np.min(np.log10(x))
#         PointLens._interpolate_1_2_x_max = np.max(np.log10(x))
#
#         with open(file_3) as file_in:
#             for line in file_in.readlines():
#                 if line[:3] == "# X":
#                     xx = np.array([float(t) for t in line.split()[2:]])
#                 if line[:3] == "# Y":
#                     yy = np.array([float(t) for t in line.split()[2:]])
#         pp = np.loadtxt(file_3)
#         PointLens._interpolate_3 = interp2d(xx, yy, pp.T, kind='cubic')
#         PointLens._interpolate_3_min_x = np.min(xx)
#         PointLens._interpolate_3_max_x = np.max(xx)
#         PointLens._interpolate_3_min_y = np.min(yy)
#         PointLens._interpolate_3_max_y = np.max(yy)
#
#         PointLens._elliptic_files_read = True
#
#     def _B_0_function(self, z):
#         """
#         calculate B_0(z) function defined in:
#
#         Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
#         https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract
#
#         Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
#         Effects from a Point-Mass Lens"
#         https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract
#
#         """
#
#         out = 4. * z / np.pi
#         def function(x): return (1. - value**2 * sin(x)**2)**.5
#
#         for (i, value) in enumerate(z):
#             if value < 1.:
#                 out[i] *= ellipe(value * value)
#             else:
#                 out[i] *= integrate.quad(function, 0.,
#                                          np.arcsin(1. / value))[0]
#         return out
#
#     def _B_1_function(self, z, B_0=None):
#         """
#         calculate B_1(z) function defined in:
#
#         Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
#         https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract
#
#         Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
#         Effects from a Point-Mass Lens"
#         https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract
#
#         """
#         if B_0 is None:
#             B_0 = self._B_0_function(z)
#
#         def function(r, theta):
#             r_2 = r * r
#             val = (1. - r_2) / (
#                 r_2 + function.arg_2 + r * function.arg_3 * cos(theta))
#             return r * sqrt(val)
#
#         def lim_0(x): return 0
#         def lim_1(x): return 1
#         rho_W_1 = 0. * z  # This equals rho * W_1().
#         for (i, zz) in enumerate(z):
#             function.arg_1 = zz
#             function.arg_2 = zz * zz
#             function.arg_3 = -2. * zz
#             rho_W_1[i] = integrate.dblquad(
#                 function, 0., 2. * np.pi, lim_0, lim_1)[0]
#
#         rho_W_1 /= np.pi
#         return B_0 - 1.5 * z * rho_W_1
#
#     def get_point_lens_finite_source_magnification(
#             self, u, pspl_magnification, direct=False):
#         """
#         Calculate magnification for point lens and finite source (for
#         a *uniform* source).  The approximation was proposed by:
#
#         `Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
#         <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_
#
#         and later the integral calculation was simplified by:
#
#         `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
#         Effects from a Point-Mass Lens"
#         <https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract>`_
#
#         This approach assumes rho is small (rho < 0.1). For larger sources
#         use :py:func:`get_point_lens_uniform_integrated_magnification`.
#
#         Parameters :
#             u: *float*, *np.array*
#                 The instantaneous source-lens separation.
#                 Multiple values can be provided.
#
#             pspl_magnification: *float*, *np.array*
#                 The point source, point lens magnification at each value of u.
#
#             direct: *boolean*
#                 Use direct calculation (slow) instead of interpolation.
#
#         Returns :
#             magnification: *float*, *np.array*
#                 The finite source source magnification.
#                 Type is the same as of u parameter.
#
#         """
#         return self._get_point_lens_finite_source_magnification(
#             u, pspl_magnification, rho=self.parameters.rho, direct=direct)
#
#     def _get_point_lens_finite_source_magnification(
#             self, u, pspl_magnification, rho, direct=False):
#         """
#         Calculate large source magnification assuming rho provided directly,
#         not as self.parameters.rho
#         """
#         z = u / rho
#         try:
#             _ = iter(z)
#         except TypeError:
#             z = np.array([z])
#
#         if direct:
#             mask = np.zeros_like(z, dtype=bool)
#         else:
#             mask = self._get_mask_B0B1_data(z)
#
#         B0 = 0. * z
#         if np.any(mask):  # Here we use interpolation.
#             B0[mask] = self._B0B1_data.interpolate_B0(z[mask])
#
#         mask = np.logical_not(mask)
#         if np.any(mask):  # Here we use direct calculation.
#             B0[mask] = self._B_0_function(z[mask])
#
#         magnification = pspl_magnification * B0
#         # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.
#         return magnification
#
#     def _get_mask_B0B1_data(self, z):
#         """
#         Get mask that desides if z is in range covered by B0B1 file
#         """
#         return self._B0B1_data.get_interpolation_mask(z)
#
#     def get_point_lens_limb_darkening_magnification(
#             self, u, pspl_magnification, gamma, direct=False):
#         """
#         calculate magnification for point lens and finite source *with
#         limb darkening*. The approximation was proposed by:
#
#         `Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
#         <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_
#
#         and later the integral calculation was simplified by:
#
#         `Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
#         Effects from a Point-Mass Lens"
#         <https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract>`_
#
#         Parameters :
#             u: *float*, *np.array*
#                 The instantaneous source-lens separation. Multiple values
#                 can be provided.
#
#             pspl_magnification: *float*, *np.array*
#                 The point source, point lens magnification at each value of u.
#
#             gamma: *float*
#                 The limb-darkening coefficient. See also
#                 :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`
#
#             direct: *boolean*
#                 Use direct calculation (very slow) instead of interpolation.
#
#         Returns :
#             magnification: *float*, *np.array*
#                 The finite source source magnification including
#                 limb-darkening. Type is the same as of u parameter.
#
#         """
#         z = u / self.parameters.rho
#         try:
#             _ = iter(z)
#         except TypeError:
#             z = np.array([z])
#
#         if direct:
#             mask = np.zeros_like(z, dtype=bool)
#         else:
#             mask = self._get_mask_B0B1_data(z)
#
#         magnification = 0. * z + pspl_magnification
#         if np.any(mask):  # Here we use interpolation.
#             B_0 = self._B0B1_data.interpolate_B0(z[mask])
#             B_0_minus_B_1 = self._B0B1_data.interpolate_B0minusB1(z[mask])
#             magnification[mask] *= (B_0 * (1. - gamma) + B_0_minus_B_1 * gamma)
#
#         mask = np.logical_not(mask)
#         if np.any(mask):  # Here we use direct calculation.
#             B_0 = self._B_0_function(z[mask])
#             B_1 = self._B_1_function(z[mask], B_0=B_0)
#             magnification[mask] *= (B_0 - gamma * B_1)
#
#         return magnification
#
#     def get_point_lens_uniform_integrated_magnification(self, u, rho):
#         """
#         Calculate magnification for the point lens and *uniform* finite source.
#         This approach works well for small and large sources
#         (e.g., rho~0.5). Uses the method presented by:
#
#         `Lee, C.-H. et al. 2009 ApJ 695, 200 "Finite-Source Effects in
#         Microlensing: A Precise, Easy to Implement, Fast, and Numerically
#         Stable Formalism"
#         <https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract>`_
#
#         Parameters :
#             u: *np.array*
#                 The instantaneous source-lens separation.
#
#             rho: *float*
#                 Source size as a fraction of the Einstein radius.
#
#         Returns :
#             magnification: *np.array*
#                 The finite source magnification.
#         """
#         n = 100
#
#         mag = np.zeros_like(u)
#
#         for i in range(len(u)):
#             if u[i] > rho:
#                 mag[i] = self._noLD_Lee09_large_u(u[i], rho, n)
#             else:
#                 mag[i] = self._noLD_Lee09_small_u(u[i], rho, n)
#         return mag
#
#     def _u_1_Lee09(self, theta, u, rho, theta_max=None):
#         """
#         Calculates Equation 4 of Lee et al. 2009.
#         The u variable is float, theta is np.ndarray.
#         """
#         if u <= rho:
#             return 0. * theta
#         out = np.zeros_like(theta)
#         mask = (theta <= theta_max)
#         if np.any(mask):
#             ucos = u * np.cos(theta[mask])
#             out[mask] = ucos - np.sqrt(rho * rho - u * u + ucos**2)
#         return out
#
#     def _u_2_Lee09(self, theta, u, rho, theta_max=None):
#         """
#         Calculates Equation 5 of Lee et al. 2009.
#         The u variable is float, theta is np.ndarray.
#         """
#         if u <= rho:
#             ucos = u * np.cos(theta)
#             return ucos + np.sqrt(rho * rho - u * u + ucos**2)
#         else:
#             out = np.zeros_like(theta)
#             mask = (theta <= theta_max)
#             if np.any(mask):
#                 ucos = u * np.cos(theta[mask])
#                 out[mask] = ucos + np.sqrt(rho * rho - u * u + ucos**2)
#             return out
#
#     def _f_Lee09(self, theta, u, rho, theta_max=None):
#         """
#         Calculates equation in text between Eq. 7 and 8 from
#         Lee et al. 2009.
#         """
#         u_1_ = self._u_1_Lee09(theta, u, rho, theta_max)
#         u_2_ = self._u_2_Lee09(theta, u, rho, theta_max)
#
#         f_u_1 = u_1_ * np.sqrt(u_1_**2 + 4.)
#         f_u_2 = u_2_ * np.sqrt(u_2_**2 + 4.)
#         return f_u_2 - f_u_1
#
#     def _noLD_Lee09_large_u(self, u, rho, n):
#         """
#         Calculates Equation 7 from Lee et al. 2009 in case u > rho.
#         """
#         if n % 2 != 0:
#             raise ValueError('internal error - odd number expected')
#         theta_max = np.arcsin(rho / u)
#         out = (u + rho) * sqrt((u + rho)**2 + 4.) - \
#             (u - rho) * sqrt((u - rho)**2 + 4.)
#         vector_1 = np.arange(1., (n / 2 - 1.) + 1)
#         vector_2 = np.arange(1., n / 2 + 1)
#         arg_1 = 2. * vector_1 * theta_max / n
#         arg_2 = (2. * vector_2 - 1.) * theta_max / n
#         out += 2. * np.sum(self._f_Lee09(arg_1, u, rho, theta_max))
#         out += 4. * np.sum(self._f_Lee09(arg_2, u, rho, theta_max))
#         out *= theta_max / (3. * np.pi * rho * rho * n)
#         return out
#
#     def _noLD_Lee09_small_u(self, u, rho, n):
#         """
#         Calculates Equation 7 from Lee et al. 2009 in case u < rho.
#         """
#         if n % 2 != 0:
#             raise ValueError('internal error - odd number expected')
#         out = (u + rho) * sqrt((u + rho)**2 + 4.) - \
#             (u - rho) * sqrt((u - rho)**2 + 4.)
#         vector_1 = np.arange(1., (n - 1.) + 1)
#         vector_2 = np.arange(1., n + 1)
#         arg_1 = vector_1 * np.pi / n
#         arg_2 = (2. * vector_2 - 1.) * np.pi / (2. * n)
#         out += 2. * np.sum(self._f_Lee09(arg_1, u, rho))
#         out += 4. * np.sum(self._f_Lee09(arg_2, u, rho))
#         out /= 2. * 3. * n * rho * rho
#         return out
#
#     def get_point_lens_LD_integrated_magnification(self, u, rho, gamma):
#         """
#         Calculate magnification for the point lens and *finite source with
#         limb-darkening*. This approach works well for small and large
#         sources (e.g., rho~0.5). Uses the method presented by:
#
#         `Lee, C.-H. et al. 2009 ApJ 695, 200 "Finite-Source Effects in
#         Microlensing: A Precise, Easy to Implement, Fast, and Numerically
#         Stable Formalism"
#         <https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract>`_
#
#         Parameters :
#             u: *np.array*
#                 The instantaneous source-lens separation.
#
#             rho: *float*
#                 Source size as a fraction of the Einstein radius.
#
#             gamma: *float*
#                 Gamma limb darkening coefficient. See also
#                 :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.
#
#         Returns :
#             magnification: *np.array*
#                 The finite source magnification.
#         """
#         n_theta = 90
#         n_u = 1000
#
#         mag = np.zeros_like(u)
#
#         for i in range(len(u)):
#             mag[i] = self._LD_Lee09(u[i], rho, gamma, n_theta, n_u)
#
#         return mag
#
#     def _LD_Lee09(self, u, rho, gamma, n_theta, n_u):
#         """
#         Calculates Equation 13 from Lee et al. 2009.
#
#         Accuracy of these calculations is on the order of 1e-4
#         """
#         theta_sub = 1.e-12
#         u_1_min = 1.e-13
#
#         if n_theta % 2 != 0:
#             raise ValueError('internal error - even number expected')
#         if n_u % 2 != 0:
#             raise ValueError('internal error - even number expected')
#
#         if u > rho:
#             theta_max = np.arcsin(rho / u)
#         else:
#             theta_max = np.pi
#         theta = np.linspace(0, theta_max - theta_sub, n_theta)
#         integrand_values = np.zeros_like(theta)
#         u_1 = self._u_1_Lee09(theta, u, rho, theta_max)
#         u_1 += u_1_min
#         u_2 = self._u_2_Lee09(theta, u, rho, theta_max)
#
#         size = (len(theta), n_u)
#         temp = np.zeros(size)
#         temp2 = (np.zeros(size).T + np.cos(theta)).T
#         for (i, (theta_, u_1_, u_2_)) in enumerate(zip(theta, u_1, u_2)):
#             temp[i] = np.linspace(u_1_, u_2_, n_u)
#         integrand = self._integrand_Lee09_v2(temp, u, temp2, rho, gamma)
#         dx = temp[:, 1] - temp[:, 0]
#         for (i, dx_) in enumerate(dx):
#             integrand_values[i] = integrate.simps(integrand[i], dx=dx_)
#         out = integrate.simps(integrand_values, dx=theta[1] - theta[0])
#         out *= 2. / (np.pi * rho**2)
#         return out
#
#     def _integrand_Lee09_v2(self, u_, u, theta_, rho, gamma):
#         """
#         Integrand in Equation 13 in Lee et al. 2009.
#
#         u_ and theta_ are np.ndarray, other parameters are scalars.
#         theta_ is in fact cos(theta_) here
#         """
#         values = 1. - (u_ * (u_ - 2. * u * theta_) + u**2) / rho**2
#         values[:, -1] = 0.
#         if values[-1, 0] < 0.:
#             values[-1, 0] = 0.
#             if values[-1, 1] < 0.:  # This sometimes happens due to rounding
#                 values[-1, 1] = .5 * values[-1, 2]  # errors above. Using
#                 # math.fsum in "values = ..." doesn't help in all cases.
#         if np.any(values < 0.):
#             if u / rho < 5.:
#                 raise ValueError(
#                     "PointLens.get_point_lens_LD_integrated_magnification() " +
#                     "unexpected error for:\nu = {:}\n".format(repr(u)) +
#                     "rho = {:}\ngamma = {:}".format(repr(rho), repr(gamma)))
#             else:
#                 message = (
#                     "PointLens.get_point_lens_LD_integrated_magnification() " +
#                     "warning! The arguments are strange: u/rho = " +
#                     "{:}.\nThere are numerical issues. You ".format(u / rho) +
#                     "can use other methods for such large u value.")
#                 warnings.warn(message, UserWarning)
#                 values[values < 0.] = 0.
#         out = 1. - gamma * (1. - 1.5 * np.sqrt(values))
#         return out * (u_**2 + 2.) / np.sqrt(u_**2 + 4.)
#
#     def get_point_lens_large_finite_source_magnification(self, u):
#         """
#         Calculate magnification for the point lens and *uniform* source.
#         This approach works well for small and large
#         sources (e.g., rho~0.5). The method was presented by:
#
#         `Witt and Mao 1994 ApJ 430, 505 "Can Lensed Stars Be Regarded as
#         Pointlike for Microlensing by MACHOs?"
#         <https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract>`_
#
#         Parameters :
#             u: *np.array*
#                 The instantaneous source-lens separation.
#
#         Returns :
#             magnification: *np.array*
#                 The finite source magnification.
#
#         """
#         out = [self._get_magnification_WM94(u_) for u_ in u]
#         return np.array(out)
#
#     def _get_magnification_WM94(self, u, rho=None):
#         """
#         Get point-lens finite-source magnification without LD.
#         """
#         if rho is None:
#             rho = self.parameters.rho
#
#         if u == rho:
#             u2 = u**2
#             a = np.pi / 2. + np.arcsin((u2 - 1.) / (u2 + 1.))
#             return (2. / u + (1. + u2) * a / u2) / np.pi
#
#         if not PointLens._elliptic_files_read:
#             self._read_elliptic_files()
#
#         a_1 = 0.5 * (u + rho) * (4. + (u - rho)**2)**.5 / rho**2
#         a_2 = -(u - rho) * (4. + 0.5 * (u**2 - rho**2))
#         a_2 /= (rho**2 * (4. + (u - rho)**2)**.5)
#         a_3 = 2. * (u - rho)**2 * (1. + rho**2)
#         a_3 /= (rho**2 * (u + rho) * (4. + (u - rho)**2)**.5)
#
#         n = 4. * u * rho / (u + rho)**2
#         k = 4. * n / (4. + (u - rho)**2)
#         # We omit sqrt, because all python packages use k^2 convention.
#
#         x_1 = self._get_ellipk(k)
#         x_2 = self._get_ellipe(k)
#         x_3 = self._get_ellip3(n, k)
#         (x_1, x_2) = (x_2, x_1)  # WM94 under Eq. 9 are inconsistent with GR80.
#
#         return (a_1 * x_1 + a_2 * x_2 + a_3 * x_3) / np.pi
#
#     def _get_ellipk(self, k):
#         """
#         Get value of elliptic integral of the first kind.
#         Use interpolation if possible.
#         """
#         x = log10(k)
#         condition_1 = (x >= PointLens._interpolate_1_2_x_min)
#         condition_2 = (x <= PointLens._interpolate_1_2_x_max)
#         if condition_1 and condition_2:
#             return PointLens._interpolate_1(x)
#         return ellipk(k)
#
#     def _get_ellipe(self, k):
#         """
#         Get value of elliptic integral of the second kind.
#         Use interpolation if possible.
#         """
#         x = log10(k)
#         condition_1 = (x >= PointLens._interpolate_1_2_x_min)
#         condition_2 = (x <= PointLens._interpolate_1_2_x_max)
#         if condition_1 and condition_2:
#             return PointLens._interpolate_2(x)
#         return ellipe(k)
#
#     def _get_ellip3(self, n, k):
#         """
#         Get value of elliptic integral of the third kind.
#         Use interpolation if possible.
#         """
#         cond_1 = (n >= PointLens._interpolate_3_min_x)
#         cond_2 = (n <= PointLens._interpolate_3_max_x)
#         cond_3 = (k >= PointLens._interpolate_3_min_y)
#         cond_4 = (k <= PointLens._interpolate_3_max_y)
#
#         if cond_1 and cond_2 and cond_3 and cond_4:
#             return PointLens._interpolate_3(n, k)[0]
#         return ellip3(n, k)
#
#     def get_point_lens_large_LD_integrated_magnification(self, u, gamma):
#         """
#         Calculate magnification for the point lens and *finite source with
#         limb-darkening*. This approach works well for small and large
#         sources (e.g., rho~0.5). Here multiple annuli
#         (each with uniform source) are used to approximate finite source with
#         limb-darkening. For uniform source calculation see:
#
#         `Witt and Mao 1994 ApJ 430, 505 "Can Lensed Stars Be Regarded as
#         Pointlike for Microlensing by MACHOs?"
#         <https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract>`_
#
#         The approximation of multiple sources in presented by, e.g.,:
#
#         `Bozza et al. 2018 MNRAS 479, 5157 "VBBINARYLENSING:
#         a public package for microlensing light-curve computation"
#         <https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract>`_
#
#         Parameters :
#             u: *np.array*
#                 The instantaneous source-lens separation.
#
#             gamma: *float*
#                 Gamma limb darkening coefficient. See also
#                 :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.
#
#         Returns :
#             magnification: *np.array*
#                 The finite source magnification.
#         """
#         n_annuli = 30  # This value could be tested better.
#
#         out = [
#             self._get_magnification_WM94_B18(u_, gamma, n_annuli) for u_ in u]
#
#         return np.array(out)
#
#     def _get_magnification_WM94_B18(self, u, gamma, n_annuli):
#         """
#         Get point-lens finite-source magnification with LD using
#         Witt & Mao 1994 approach and equations 16-19 from Bozza et al. 2018.
#         """
#         n_annuli += 1  # It's easier to have r=0 ring as well.
#         annuli = np.linspace(0, 1., n_annuli)
#         r2 = annuli**2
#
#         magnification = np.zeros(n_annuli)
#         for (i, a) in enumerate(annuli):
#             if i == 0:
#                 continue
#             magnification[i] = self._get_magnification_WM94(
#                 u=u, rho=a * self.parameters.rho)
#
#         cumulative_profile = gamma + (1. - gamma) * r2 - gamma * (1. - r2)**1.5
#         d_cumulative_profile = cumulative_profile[1:] - cumulative_profile[:-1]
#         d_r2 = r2[1:] - r2[:-1]
#         temp = magnification * r2
#         d_mag_r2 = temp[1:] - temp[:-1]
#         out = np.sum(d_mag_r2 * d_cumulative_profile / d_r2)
#         return out
#

class PointSourcePointLensMagnification():
    """
    Equations for calculating point-source--point-lens magnification and
    its derivatives.

    Keywords :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`

    """

    def __init__(self, trajectory=None):
        if not isinstance(trajectory, mm.Trajectory):
            raise TypeError(
                "PointLens* argument has to be of Trajectory type, not " +
                str(type(trajectory)))

        self.trajectory = trajectory

        self._u_2 = None # u^2
        self._u = None

        self._pspl_magnification = None
        self._magnification = None

    def get_pspl_magnification(self):
        """
        This is Paczynski equation, i.e., point-source--point-lens (PSPL)
        magnification.
        Arguments :
         Parameters :
            u: *np.array*
                The instantaneous source-lens separation.
        Returns :
            pspl_magnification: *float* or *np.ndarray*
                The point-source--point-lens magnification for each point
                specified by `u`.
        """
        self._pspl_magnification = \
            (self.u_2 + 2.) / np.sqrt(self.u_2 * (self.u_2 + 4.))

        return self._pspl_magnification

    def get_magnification(self):
        """
        Arguments :
         Parameters :
            u: *np.array*
                The instantaneous source-lens separation.
        Returns :
            magnification: *float* or *np.ndarray*
                The magnification for each point
                specified by `u.
        """
        self._magnification = self.get_pspl_magnification()
        return self._magnification

    def get_d_A_d_params(self, parameters):
        """
        Calculate d A / d parameters for a point lens model.

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            dA_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """
        d_A_d_params = {}
        d_u_d_params = self.get_d_u_d_params(parameters)
        d_A_d_u = self.get_d_A_d_u()
        for key in parameters:
            d_A_d_params[key] = d_A_d_u * d_u_d_params[key]

        return d_A_d_params

    def get_d_u_d_params(self, parameters):
        """
        Calculate d u / d parameters

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            du_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """
        # Setup
        d_u_d_params = {param: 0 for param in parameters}
        as_dict = self.trajectory.parameters.as_dict()

        # Calculate derivatives
        d_u_d_x = self.trajectory.x / self.u_
        d_u_d_y = self.trajectory.y / self.u_
        dt = self.trajectory.times - as_dict['t_0']

        # Exactly 2 out of (u_0, t_E, t_eff) must be defined and
        # gradient depends on which ones are defined.
        t_E = self.trajectory.parameters.t_E
        t_eff = self.trajectory.parameters.t_eff
        if 't_eff' not in as_dict:
            d_u_d_params['t_0'] = -d_u_d_x / t_E
            d_u_d_params['u_0'] = d_u_d_y
            d_u_d_params['t_E'] = d_u_d_x * -dt / t_E**2
        elif 't_E' not in as_dict:
            d_u_d_params['t_0'] = -d_u_d_x * as_dict['u_0'] / t_eff
            d_u_d_params['u_0'] = (d_u_d_y + d_u_d_x * dt / t_eff)
            d_u_d_params['t_eff'] = (d_u_d_x * -dt * as_dict['u_0'] / t_eff**2)
        elif 'u_0' not in as_dict:
            d_u_d_params['t_0'] = -d_u_d_x / t_E
            d_u_d_params['t_E'] = (d_u_d_x * dt - d_u_d_y * t_eff) / t_E**2
            d_u_d_params['t_eff'] = d_u_d_y / t_E
        else:
            raise KeyError(
                'Something is wrong with ModelParameters in ' +
                'FitData.calculate_chi2_gradient():\n', as_dict)

        # Below we deal with parallax only.
        if 'pi_E_N' in parameters or 'pi_E_E' in parameters:
            delta_N = self.trajectory.parallax_delta_N_E['N']
            delta_E = self.trajectory.parallax_delta_N_E['E']

            d_u_d_params['pi_E_N'] = d_u_d_x * delta_N + d_u_d_y * delta_E
            d_u_d_params['pi_E_E'] = d_u_d_x * delta_E - d_u_d_y * delta_N

        return d_u_d_params

    def get_d_A_d_u(self):
        """
        Calculate dA/du for PSPL point-source--point-lens model.

        No parameters.

        Returns :
            dA_du: *np.ndarray*
                Derivative dA/du evaluated at each epoch.
        """
        d_A_d_u = -8. / (self.u_2 * (self.u_2 + 4) * np.sqrt(self.u_2 + 4))
        return d_A_d_u

    @property
    def pspl_magnification(self):
        """
        *np.ndarray*

        Point-source--point-lens magnification for each epoch.
        """
        if self._pspl_magnification is None:
            self.get_pspl_magnification()

        return self._pspl_magnification

    @property
    def magnification(self):
        """
        *np.ndarray*

        Magnification for each epoch.
        """
        return self._magnification

    @property
    def u_(self):
        """
        *np.ndarray*

        Magnitude of lens-source separation for each epoch.
        """
        if self._u is None:
            self._u = np.sqrt(self.u_2)

        return self._u

    @property
    def u_2(self):
        """
        *np.ndarray*

        Square of the magnitude of lens-source separation for each epoch.
        """
        if self._u_2 is None:
            self._u_2 = self.trajectory.x ** 2 + self.trajectory.y ** 2

        return self._u_2


class FiniteSourceUniformGould94Magnification(
    PointSourcePointLensMagnification):
    """
    Equations for calculating finite-source--point-lens magnification and
    its derivatives following the `Gould 1994 ApJ, 421L, 71
    <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_
    prescription assuming a *uniform* (and circular) source.

    Keywords :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`

    """

    def __init__(self, direct=False, **kwargs):
        PointSourcePointLensMagnification.__init__(self, **kwargs)

        self.direct = direct
        self._B0B1_data = mm.B0B1Utils()

        self._z = None
        self._b0 = None
        self._db0 = None

    def get_magnification(self):
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
        pspl_magnification = self.get_pspl_magnification()
        self._magnification = pspl_magnification * self.b0
        # More accurate calculations can be performed - see Yoo+04 eq. 11 & 12.

        return self._magnification

    def _B_0_function(self, mask=None):
        """
        calculate B_0(z) function defined in:
        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract
        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract
        """
        if mask is not None:
            z = self.z_[mask]
        else:
            z = self.z_

        out = 4. * z / np.pi

        def function(x):
            return (1. - value ** 2 * sin(x) ** 2) ** .5

        for (i, value) in enumerate(z):
            if value < 1.:
                out[i] *= ellipe(value * value)
            else:
                out[i] *= integrate.quad(function, 0.,
                                         np.arcsin(1. / value))[0]
        return out

    def _get_fspl_deriv_factor(self):
        factor = self.pspl_magnification * self.db0
        factor /= self.trajectory.parameters.rho
        factor += self.get_d_A_d_u() * self.b0

        return factor

    def get_d_A_d_params(self, parameters):
        """
        Return the gradient of the magnification with respect to the
        FSPL parameters for each epoch.

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns:
            du_dparams: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """

        d_u_d_params = PointSourcePointLensMagnification.get_d_u_d_params(
            self, parameters)

        factor = self._get_fspl_deriv_factor()

        d_A_d_params = {}
        for key in parameters:
            if key == 'rho':
                d_A_d_params[key] = self.get_d_A_d_rho()
            else:
                d_A_d_params[key] = d_u_d_params[key] * factor

        return d_A_d_params

    def get_d_A_d_rho(self):
        """
        Return the derivative of the magnification with respect to rho
        for each epoch.

        No parameters.

        Returns :
            dA_drho: *np.ndarray*
                Derivative dA/drho evaluated at each epoch.

        """
        d_A_d_rho = np.ones(len(self.trajectory.times))
        d_A_d_rho *= self.pspl_magnification
        d_A_d_rho *= -self.u_ / self.trajectory.parameters.rho**2
        d_A_d_rho *= self.db0

        return d_A_d_rho

    @property
    def z_(self):
        """ Magnitude of lens-source separation scaled to rho for each epoch."""
        if self._z is None:
            self._z = self.u_ / self.trajectory.parameters.rho

        return self._z

    @property
    def b0(self):
        """
        *np.ndarray*

        Return the value of B_0(z) function for each epoch.
        """
        if self._b0 is None:
            if self.direct:
                mask = np.zeros_like(self.z_, dtype=bool)
            else:
                mask = self._B0B1_data.get_interpolation_mask(self.z_)

            self._b0 = 0. * self.z_
            if np.any(mask):  # Here we use interpolation.
                self._b0[mask] = self._B0B1_data.interpolate_B0(self.z_[mask])

            mask = np.logical_not(mask)
            if np.any(mask):  # Here we use direct calculation.
                self._b0[mask] = self._B_0_function(mask)

        return self._b0

    @property
    def db0(self):
        """
        *np.ndarray*

        Retrieve derivative of B_0(z) function for each epoch.
        """
        if self._db0 is None:
            if self.direct:
                raise NotImplementedError(
                    'B0 derivatives not implemented for direct method.')
            else:
                mask = self._B0B1_data.get_interpolation_mask(self.z_)

            self._db0 = 0. * self.z_
            if np.any(mask):  # Here we use interpolation.
                self._db0[mask] = self._B0B1_data.interpolate_B0prime(
                    self.z_[mask])

        return self._db0


class FiniteSourceLDYoo04Magnification(FiniteSourceUniformGould94Magnification):

    def __init__(self, gamma=None, **kwargs):
        FiniteSourceUniformGould94Magnification.__init__(self, **kwargs)

        self._gamma = gamma
        self._b1 = None
        self._db1 = None

    def get_magnification(self):
        FiniteSourceUniformGould94Magnification.get_magnification(self)
        self._magnification -= self.pspl_magnification * self.b1 * self._gamma

        return self._magnification

    def _get_fspl_deriv_factor(self):
        factor = self.pspl_magnification * (self.db0 - self._gamma * self.db1)
        factor /= self.trajectory.parameters.rho
        factor += self.get_d_A_d_u() * (self.b0 - self._gamma * self.b1)

        return factor

    def get_d_A_d_rho(self):
        """
        Return the derivative of the magnification with respect to rho
        for each epoch.

        No parameters.

        Returns :
            dA_drho: *np.ndarray*
                Derivative dA/drho evaluated at each epoch.

        """
        d_A_d_rho = np.ones(len(self.trajectory.times))
        d_A_d_rho *= self.pspl_magnification
        d_A_d_rho *= -self.u_ / self.trajectory.parameters.rho**2
        d_A_d_rho *= (self.db0 - self._gamma * self.db1)
        return d_A_d_rho

    def _B_1_function(self, mask=None):
        """
        calculate B_1(z) function defined in:

        Gould A. 1994 ApJ 421L, 71 "Proper motions of MACHOs"
        https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract

        Yoo J. et al. 2004 ApJ 603, 139 "OGLE-2003-BLG-262: Finite-Source
        Effects from a Point-Mass Lens"
        https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract

        """
        if mask is not None:
            z = self.z_[mask]
        else:
            z = self.z_

        def function(r, theta):
            r_2 = r * r
            val = (1. - r_2) / (
                r_2 + function.arg_2 + r * function.arg_3 * cos(theta))
            return r * sqrt(val)

        def lim_0(x): return 0
        def lim_1(x): return 1

        rho_W_1 = 0. * z  # This equals rho * W_1().
        for (i, zz) in enumerate(z):
            function.arg_1 = zz
            function.arg_2 = zz * zz
            function.arg_3 = -2. * zz
            rho_W_1[i] = integrate.dblquad(
                function, 0., 2. * np.pi, lim_0, lim_1)[0]

        rho_W_1 /= np.pi

        return self.b0 - 1.5 * z * rho_W_1

    @property
    def b1(self):
        """
        *np.ndarray*

        Return the value of B_1(z) function for each epoch.
        """
        if self._b1 is None:
            if self.direct:
                mask = np.zeros_like(self.z_, dtype=bool)
            else:
                mask = self._B0B1_data.get_interpolation_mask(self.z_)

            self._b1 = 0. * self.z_
            if np.any(mask):  # Here we use interpolation.
                self._b1[mask] = self._B0B1_data.interpolate_B1(self.z_[mask])

            mask = np.logical_not(mask)
            if np.any(mask):  # Here we use direct calculation.
                self._b1[mask] = self._B_1_function(mask)

        return self._b1

    @property
    def db1(self):
        """
        *np.ndarray*

        Retrieve derivative of B_1(z) function for each epoch.
        """
        if self._db1 is None:
            if self.direct:
                raise NotImplementedError(
                    'B0 derivatives not implemented for direct method.')
            else:
                mask = self._B0B1_data.get_interpolation_mask(self.z_)

            self._db1 = 0. * self.z_
            if np.any(mask):  # Here we use interpolation.
                self._db1[mask] = self._B0B1_data.interpolate_B1prime(
                    self.z_[mask])


        return self._db1


class FiniteSourceUniformWittMao94Magnification(
    PointSourcePointLensMagnification):
    """
    Calculate magnification for the point lens and *uniform* source.
    This approach works well for small and large
    sources (e.g., rho~0.5). The method was presented by:

    `Witt and Mao 1994 ApJ 430, 505 "Can Lensed Stars Be Regarded as
    Pointlike for Microlensing by MACHOs?"
    <https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract>`_
    """

    def __init__(self, **kwargs):
        PointSourcePointLensMagnification.__init__(self, **kwargs)

        self._B0B1_data = mm.B0B1Utils()

    def get_magnification(self):
        out = [self._get_magnification_WM94(u_) for u_ in self.u_]
        return np.array(out)

    def _get_magnification_WM94(self, u):
        """
        Get point-lens finite-source magnification without LD.
        """
        rho = self.trajectory.parameters.rho

        if u == rho:
            u2 = u**2
            a = np.pi / 2. + np.arcsin((u2 - 1.) / (u2 + 1.))
            return (2. / u + (1. + u2) * a / u2) / np.pi

        a_1 = 0.5 * (u + rho) * (4. + (u - rho)**2)**.5 / rho**2
        a_2 = -(u - rho) * (4. + 0.5 * (u**2 - rho**2))
        a_2 /= (rho**2 * (4. + (u - rho)**2)**.5)
        a_3 = 2. * (u - rho)**2 * (1. + rho**2)
        a_3 /= (rho**2 * (u + rho) * (4. + (u - rho)**2)**.5)

        n = 4. * u * rho / (u + rho)**2
        k = 4. * n / (4. + (u - rho)**2)
        # We omit sqrt, because all python packages use k^2 convention.

        x_1 = self._get_ellipk(k)
        x_2 = self._get_ellipe(k)
        x_3 = self._get_ellip3(n, k)
        (x_1, x_2) = (x_2, x_1)  # WM94 under Eq. 9 are inconsistent with GR80.

        return (a_1 * x_1 + a_2 * x_2 + a_3 * x_3) / np.pi

    def _get_ellipk(self, k):
        """
        Get value of elliptic integral of the first kind.
        Use interpolation if possible.
        """
        x = log10(k)
        condition_1 = (x >= mm.EllipUtils._interpolate_1_2_x_min)
        condition_2 = (x <= mm.EllipUtils._interpolate_1_2_x_max)
        if condition_1 and condition_2:
            return mm.EllipUtils._interpolate_1(x)

        return ellipk(k)

    def _get_ellipe(self, k):
        """
        Get value of elliptic integral of the second kind.
        Use interpolation if possible.
        """
        x = log10(k)
        condition_1 = (x >= mm.EllipUtils._interpolate_1_2_x_min)
        condition_2 = (x <= mm.EllipUtils._interpolate_1_2_x_max)
        if condition_1 and condition_2:
            return mm.EllipUtils._interpolate_2(x)

        return ellipe(k)

    def _get_ellip3(self, n, k):
        """
        Get value of elliptic integral of the third kind.
        Use interpolation if possible.
        """
        cond_1 = (n >= mm.EllipUtils._interpolate_3_min_x)
        cond_2 = (n <= mm.EllipUtils._interpolate_3_max_x)
        cond_3 = (k >= mm.EllipUtils._interpolate_3_min_y)
        cond_4 = (k <= mm.EllipUtils._interpolate_3_max_y)

        if cond_1 and cond_2 and cond_3 and cond_4:
            return mm.EllipUtils._interpolate_3(n, k)[0]

        return ellip3(n, k)

    def get_d_A_d_params(self, parameters):
        """
        Calculate d A / d parameters for a point lens model.

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            dA_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """
        raise NotImplementedError(
            'Derivative calculations Not Implemented for WittMao94')

    def get_d_u_d_params(self, parameters):
        """
        Calculate d u / d parameters

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            du_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """
        raise NotImplementedError(
            'Derivative calculations Not Implemented for WittMao94')

    def get_d_A_d_u(self):
        """
        Calculate dA/du for PSPL point-source--point-lens model.

        No parameters.

        Returns :
            dA_du: *np.ndarray*
                Derivative dA/du evaluated at each epoch.
        """
        raise NotImplementedError(
            'Derivative calculations Not Implemented for WittMao94')

    def get_d_A_d_rho(self):
        """
        Return the derivative of the magnification with respect to rho
        for each epoch.

        No parameters.

        Returns :
            dA_drho: *np.ndarray*
                Derivative dA/drho evaluated at each epoch.

        """
        raise NotImplementedError(
            'Derivative calculations Not Implemented for WittMao94')

class FiniteSourceLDWittMao94Magnification(
    FiniteSourceUniformWittMao94Magnification):
    """
    Calculate magnification for the point lens and *finite source with
    limb-darkening*. This approach works well for small and large
    sources (e.g., rho~0.5). Here multiple annuli
    (each with uniform source) are used to approximate finite source with
    limb-darkening. For uniform source calculation see:

    `Witt and Mao 1994 ApJ 430, 505 "Can Lensed Stars Be Regarded as
    Pointlike for Microlensing by MACHOs?"
    <https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract>`_

    The approximation of multiple sources in presented by, e.g.,:

    `Bozza et al. 2018 MNRAS 479, 5157 "VBBINARYLENSING:
    a public package for microlensing light-curve computation"
    <https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract>`_

    Parameters :
        u: *np.array*
            The instantaneous source-lens separation.

        gamma: *float*
            Gamma limb darkening coefficient. See also
            :py:class:`~MulensModel.limbdarkeningcoeffs.LimbDarkeningCoeffs`.

    Returns :
        magnification: *np.array*
            The finite source magnification.
    """

    def __init__(self, gamma=None, **kwargs):
        FiniteSourceUniformWittMao94Magnification.__init__(self, **kwargs)

        self._gamma = gamma
        self.n_annuli = 30  # This value could be tested better.

    def get_magnification(self):
        out = [
            self._get_magnification_WM94_B18(u_)
            for u_ in self.u_]

        return np.array(out)

    def _get_magnification_WM94_B18(self, u):
        """
        Get point-lens finite-source magnification with LD using
        Witt & Mao 1994 approach and equations 16-19 from Bozza et al. 2018.
        """
        n_annuli = self.n_annuli + 1  # It's easier to have r=0 ring as well.
        annuli = np.linspace(0, 1., n_annuli)
        r2 = annuli**2

        magnification = np.zeros(n_annuli)
        for (i, a) in enumerate(annuli):
            if i == 0:
                continue
            magnification[i] = self._get_magnification_WM94(
                u=u, rho=a * self.trajectory.parameters.rho)

        cumulative_profile = (self.gamma + (1. - self.gamma) * r2 -
                              self.gamma * (1. - r2)**1.5)
        d_cumulative_profile = cumulative_profile[1:] - cumulative_profile[:-1]
        d_r2 = r2[1:] - r2[:-1]
        temp = magnification * r2
        d_mag_r2 = temp[1:] - temp[:-1]
        out = np.sum(d_mag_r2 * d_cumulative_profile / d_r2)

        return out


class FiniteSourceUniformLee09Magnification(PointSourcePointLensMagnification):
    pass


class FiniteSourceLDLee09Magnification(PointSourcePointLensMagnification):
    pass
