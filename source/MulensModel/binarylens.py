import sys
import os
import ctypes
import glob
import warnings
import numpy as np
from math import fsum, sqrt

import MulensModel as mm
try:
    import MulensModel.VBBL as mm_vbbl
except Exception:
    _vbbl_wrapped = False
else:
    _vbbl_wrapped = True
try:
    import MulensModel.AdaptiveContouring as mm_ac
except Exception:
    _adaptive_contouring_wrapped = False
else:
    _adaptive_contouring_wrapped = True


def _try_load(path, name):
    """
    Try loading compiled C library.
    Input is *str* or *list* of *str*.
    """
    if isinstance(path, str):
        path = [path]
    for path_ in path:
        try:
            out = ctypes.cdll.LoadLibrary(path_)
        except OSError:
            print("WARNING - File not loaded:", path_)
            print("Everything should work except:", name)
            pass
        else:
            return out
    return None


def _get_path_2(name_1, name_2):
    """convenience function"""
    module_path = os.path.abspath(__file__)
    for i in range(3):
        module_path = os.path.dirname(module_path)
    return os.path.join(module_path, 'source', name_1, name_2)


def _import_compiled_VBBL():
    """try importing manually compiled VBBL package"""
    vbbl = _try_load(
        _get_path_2('VBBL', "VBBinaryLensingLibrary_wrapper.so"), "VBBL")
    _vbbl_wrapped = (vbbl is not None)
    if not _vbbl_wrapped:
        return (_vbbl_wrapped, None, None, None, None)

    vbbl.VBBinaryLensing_BinaryMagDark.argtypes = 7 * [ctypes.c_double]
    vbbl.VBBinaryLensing_BinaryMagDark.restype = ctypes.c_double

    vbbl.VBBinaryLensing_BinaryMag0.argtypes = 7 * [ctypes.c_double]
    vbbl.VBBinaryLensing_BinaryMag0.restype = ctypes.c_double

    vbbl.VBBL_SG12_5.argtypes = 12 * [ctypes.c_double]
    vbbl.VBBL_SG12_5.restype = np.ctypeslib.ndpointer(
        dtype=ctypes.c_double, shape=(10,))

    vbbl.VBBL_SG12_9.argtypes = 20 * [ctypes.c_double]
    vbbl.VBBL_SG12_9.restype = np.ctypeslib.ndpointer(
            dtype=ctypes.c_double, shape=(18,))

    return (_vbbl_wrapped,
            vbbl.VBBinaryLensing_BinaryMagDark, vbbl.VBBL_SG12_5, 
            vbbl.VBBinaryLensing_BinaryMag0, vbbl.VBBL_SG12_9)


def _import_compiled_AdaptiveContouring():
    """try importing manually compiled AdaptiveContouring package"""
    ac = "AdaptiveContouring"
    adaptive_contour = _try_load(_get_path_2(ac, ac + "_wrapper.so"), ac)
    _adaptive_contouring_wrapped = (adaptive_contour is not None)
    if not _adaptive_contouring_wrapped:
        return (_adaptive_contouring_wrapped, None)
    adaptive_contour.Adaptive_Contouring_Linear.argtypes = (
        8 * [ctypes.c_double])
    adaptive_contour.Adaptive_Contouring_Linear.restype = ctypes.c_double
    return (_adaptive_contouring_wrapped,
            adaptive_contour.Adaptive_Contouring_Linear)



# Check import and try manually compiled versions.
if _vbbl_wrapped:
    _vbbl_binary_mag_dark = mm_vbbl.VBBinaryLensing_BinaryMagDark
    _vbbl_binary_mag_0 = mm_vbbl.VBBinaryLensing_BinaryMag0
    _vbbl_SG12_5 = mm_vbbl.VBBL_SG12_5
    _vbbl_SG12_9 = mm_vbbl.VBBL_SG12_9
else:
    out = _import_compiled_VBBL()
    _vbbl_wrapped = out[0]
    _vbbl_binary_mag_dark = out[1]
    _vbbl_SG12_5 = out[2]
    _vbbl_binary_mag_0 = out[3]
    _vbbl_SG12_9 = out[4]

if not _vbbl_wrapped:
    _solver = 'numpy'
else:
    _solver = 'Skowron_and_Gould_12'
if _adaptive_contouring_wrapped:
    _adaptive_contouring_linear = mm_ac.Adaptive_Contouring_Linear
else:
    out = _import_compiled_AdaptiveContouring()
    _adaptive_contouring_wrapped = out[0]
    _adaptive_contouring_linear = out[1]

class BinaryLens(object):
    """
    The binary lens equation - its solutions, images, parities,
    magnifications, etc.

    The binary lens equation is a 5th order complex polynomial.

    Attributes :
        mass_1: *float*
            mass of the primary (left-hand object) as a fraction of
            the total mass.

        mass_2: *float*
            mass of the secondary (right-hand object) as a fraction of the
            total mass.

        separation: *float*
            separation between the two bodies as a fraction of the Einstein
            ring.

    Note: mass_1 and mass_2 may be defined as a fraction of some other
    mass than the total mass. This is possible but not recommended -
    make sure you know what you're doing before you start using this
    possibility.

    """
    def __init__(self, mass_1=None, mass_2=None, separation=None):
        self.mass_1 = float(mass_1)  # This speeds-up code for np.float input.
        self.mass_2 = float(mass_2)
        self.separation = float(separation)
        self._total_mass = None
        self._mass_difference = None
        self._position_z1 = None
        self._position_z2 = None
        self._last_polynomial_input = None
        self._solver = _solver
        self._use_planet_frame = True

    def _calculate_variables(self, source_x, source_y):
        """calculates values of constants needed for polynomial coefficients"""
        self._total_mass = 0.5 * (self.mass_1 + self.mass_2)
        # This is total_mass in WM95 paper.

        self._mass_difference = 0.5 * (self.mass_2 - self.mass_1)
        self._zeta = source_x + source_y * 1.j
        if self._use_planet_frame:
            self._position_z1 = -self.separation + 0.j
            self._position_z2 = 0. + 0.j
        else:
            self._position_z1 = -0.5 * self.separation + 0.j
            self._position_z2 = 0.5 * self.separation + 0.j

    def _get_polynomial(self, source_x, source_y):
        """get polynomial coefficients"""
        if self._use_planet_frame:
            return self._get_polynomial_planet_frame(source_x, source_y)
        else:
            return self._get_polynomial_WM95(source_x, source_y)

    def _get_polynomial_WM95(self, source_x, source_y):
        """
        calculate coefficients of the polynomial in geometric center frame
        """
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff

        pos_z1 = self._position_z1

        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2

        zeta = self._zeta
        zeta_conj = zeta.conjugate()
        zeta_conj_pow2 = zeta_conj * zeta_conj

        # Calculate the coefficients of the 5th order complex polynomial
        coeff_5 = mm.Utils.complex_fsum([z1_pow2, -zeta_conj_pow2])
        coeff_4 = mm.Utils.complex_fsum(
            [-2. * total_m * zeta_conj,
             zeta * zeta_conj_pow2, -2. * m_diff * pos_z1,
             -zeta * z1_pow2])
        coeff_3 = mm.Utils.complex_fsum(
            [4. * total_m * zeta * zeta_conj,
             4. * m_diff * zeta_conj * pos_z1,
             2. * zeta_conj_pow2 * z1_pow2, -2. * z1_pow4])
        coeff_2 = mm.Utils.complex_fsum(
            [4. * total_m_pow2 * zeta,
             4. * total_m * m_diff * pos_z1,
             -4. * m_diff * zeta * zeta_conj * pos_z1,
             -2. * zeta * zeta_conj_pow2 * z1_pow2,
             4. * m_diff * z1_pow3, 2. * zeta * z1_pow4])
        coeff_1 = mm.Utils.complex_fsum(
            [-8. * total_m * m_diff * zeta * pos_z1,
             -4. * m_diff_pow2 * z1_pow2,
             -4. * total_m_pow2 * z1_pow2,
             -4. * total_m * zeta * zeta_conj * z1_pow2,
             -4. * m_diff * zeta_conj * z1_pow3,
             -zeta_conj_pow2 * z1_pow4, z1_pow3 * z1_pow3])
        coeff_0 = mm.Utils.complex_fsum(
            [4. * m_diff_pow2 * zeta,
             4. * total_m * m_diff * pos_z1,
             4. * m_diff * zeta * zeta_conj * pos_z1,
             2. * total_m * zeta_conj * z1_pow2,
             zeta * zeta_conj_pow2 * z1_pow2,
             -2. * m_diff * z1_pow3 - zeta * z1_pow4])
        coeff_0 *= z1_pow2

        # Return the coefficients of the polynomial
        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5]
        return np.array(coeffs_list).reshape(6)

    def _get_polynomial_planet_frame(self, source_x, source_y):
        """calculate coefficients of the polynomial in planet frame"""
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass

        m_diff = self._mass_difference

        zeta = self._zeta
        zeta_conj = zeta.conjugate()

        c_sum = mm.Utils.complex_fsum

        z1 = self._position_z1

        coeff_5 = c_sum([z1, -zeta_conj]) * zeta_conj
        coeff_4 = c_sum([
            (-m_diff + total_m) * z1,
            -c_sum([2. * total_m, z1 * c_sum([2. * z1, zeta])]) * zeta_conj,
            c_sum([2. * z1 + zeta]) * zeta_conj**2
            ])
        coeff_3 = c_sum([
            z1 * c_sum([m_diff * z1, -total_m * c_sum([z1, 2. * zeta])]),
            zeta_conj * c_sum([
                2. * m_diff * z1,
                c_sum([2. * total_m, z1**2]) * c_sum([z1, 2. * zeta])
                ]),
            -z1 * c_sum([z1, 2. * zeta]) * zeta_conj**2
            ])
        coeff_2 = c_sum([
            m_diff * z1 * c_sum([2. * total_m, z1 * zeta]),
            total_m * c_sum([
                -2. * total_m * z1, 4. * total_m * zeta, 3. * z1**2 * zeta]),
            -z1 * zeta_conj * c_sum([
                zeta * c_sum([6. * total_m, z1**2]),
                2. * m_diff * c_sum([z1, zeta])
                ]),
            z1**2 * zeta * zeta_conj**2
            ])
        coeff_1 = -z1 * (m_diff + total_m) * c_sum([
            m_diff * z1, -total_m * z1, 4. * total_m * zeta, z1**2 * zeta,
            -2. * z1 * zeta * zeta_conj
            ])
        coeff_0 = (m_diff + total_m)**2 * z1**2 * zeta

        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5]
        return np.array(coeffs_list).reshape(6)

    def _get_polynomial_roots(self, source_x, source_y):
        """roots of the polynomial"""
        polynomial_input = [self.mass_1, self.mass_2, self.separation,
                            source_x, source_y]

        if polynomial_input == self._last_polynomial_input:
            return self._polynomial_roots

        polynomial = self._get_polynomial(
            source_x=source_x, source_y=source_y)

        np_polyroots = np.polynomial.polynomial.polyroots
        if self._solver == 'numpy':
            self._polynomial_roots = np_polyroots(polynomial)
        elif self._solver == 'Skowron_and_Gould_12':
            args = polynomial.real.tolist() + polynomial.imag.tolist()
            try:
                out = _vbbl_SG12_5(*args)
            except ValueError as err:
                err2 = "\n\nSwitching from Skowron & Gould 2012 to numpy"
                warnings.warn(str(err) + err2, UserWarning)
                self._solver = 'numpy'
                self._polynomial_roots = np_polyroots(polynomial)
            else:
                self._polynomial_roots = np.array([
                    out[0]+out[5]*1.j, out[1]+out[6]*1.j, out[2]+out[7]*1.j,
                    out[3]+out[8]*1.j, out[4]+out[9]*1.j])
        else:
            raise ValueError('Unknown solver: {:}'.format(self._solver))
        self._last_polynomial_input = polynomial_input

        return self._polynomial_roots

    def _polynomial_roots_ok(
            self, source_x, source_y, return_distances=False):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots(
            source_x=source_x, source_y=source_y)

        # Two lines below are simplified assuming
        # self._position_z1.imag = 0 and same for z2.
        roots_conj = np.conjugate(roots)
        solutions = (self._zeta +
                     self.mass_1 / (roots_conj - self._position_z1) +
                     self.mass_2 / (roots_conj - self._position_z2))
        # This backs-up the lens equation.

        out = []
        distances = []
        for (i, root) in enumerate(roots):
            distances_from_root = abs((solutions-root)**2)
            min_distance_arg = np.argmin(distances_from_root)

            if i == min_distance_arg:
                out.append(root)
                distances.append(distances_from_root[min_distance_arg])
            # The values in distances[] are a diagnostic on how good the
            # numerical accuracy is.

        # If the lens equation is solved correctly, there should be
        # either 3 or 5 solutions (corresponding to 3 or 5 images)
        if len(out) not in [3, 5]:
            msg = ("Wrong number of solutions to the lens equation of binary" +
                   " lens.\nGot {:} and expected 3 or 5.\nThe parameters " +
                   "(m1, m2, s, source_x, source_y, solver) are:\n" +
                   "{:} {:} {:} {:} {:}  {:}\n\n" +
                   "Consider using 'point_source_point_lens' method for " +
                   "epochs when the source is very far from the lens. Note " +
                   "that it's different from 'point_source' method.")
            txt = msg.format(
                len(out), repr(self.mass_1), repr(self.mass_2),
                repr(self.separation), repr(source_x), repr(source_y),
                self._solver)

            if self._solver != "Skowron_and_Gould_12":
                txt += (
                    "\n\nYou should switch to using Skowron_and_Gould_12" +
                    " polynomial root solver. It is much more accurate than " +
                    "numpy.polynomial.polynomial.polyroots(). " +
                    "Skowron_and_Gould_12 method is selected in automated " +
                    "way if VBBL is imported properly.")
            distance = sqrt(source_x**2 + source_y**2)
            if (self.mass_2 > 1.e-6 * self.mass_1 and
                    (distance < 15. or distance < 2. * self.separation)):
                txt += ("\n\nThis is surprising error - please contact code " +
                        "authors and provide the above error message.")
            elif distance > 200.:
                txt += ("\n\nYou try to calculate magnification at huge " +
                        "distance from the source and this is causing an " +
                        "error.")
            txt += "\nMulensModel version: {:}".format(mm.__version__)

            raise ValueError(txt)

        if return_distances:
            return (np.array(out), np.array(distances))
        else:
            return np.array(out)

    def _jacobian_determinant_ok(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._polynomial_roots_ok(
                                   source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        add_1 = self.mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self.mass_2 / (self._position_z2 - roots_ok_bar)**2
        derivative = add_1 + add_2

        return 1. - derivative * np.conjugate(derivative)

    def _signed_magnification(self, source_x, source_y):
        """signed magnification for each image separately"""
        return 1. / self._jacobian_determinant_ok(
                source_x=source_x, source_y=source_y)

    def _point_source(self, source_x, source_y):
        """calculate point source magnification"""
        signed_magnification = self._signed_magnification(
            source_x=source_x, source_y=source_y)
        return fsum(abs(signed_magnification))

    def _point_source_Witt_Mao_95(self, source_x, source_y):
        """calculate point source magnification"""
        return self._point_source(source_x=source_x, source_y=source_y)

    def point_source_magnification(self, source_x, source_y):
        """
        Calculate point source magnification for given position. The
        origin of the coordinate system is at the center of mass and
        both masses are on X axis with higher mass at negative X; this
        means that the higher mass is at (X, Y)=(-s*q/(1+q), 0) and
        the lower mass is at (s/(1+q), 0).

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

        Returns :
            magnification: *float*
                Point source magnification.
        """
        if self._use_planet_frame:
            x_shift = -self.mass_1 / (self.mass_1 + self.mass_2)
        else:
            x_shift = self.mass_2 / (self.mass_1 + self.mass_2) - 0.5
        x_shift *= self.separation
        # We need to add this because in order to shift to correct frame.
        return self._point_source_Witt_Mao_95(
                source_x=float(source_x)+x_shift, source_y=float(source_y))
        # Casting to float speeds-up code for np.float input.

    def _get_magnification_w_plus(self, source_x, source_y, radius,
                                  magnification_center=None):
        """Evaluates Gould (2008) eq. 7"""
        dx = [1., 0., -1., 0.]
        dy = [0., 1., 0., -1.]
        out = []
        for (i, dxval) in enumerate(dx):
            x = source_x + dxval * radius
            y = source_y + dy[i] * radius
            out.append(self.point_source_magnification(
                                              source_x=x, source_y=y))
        if magnification_center is None:
            magnification_center = self.point_source_magnification(
                                    source_x=source_x, source_y=source_y)
        return 0.25 * fsum(out) - magnification_center

    def _get_magnification_w_times(self, source_x, source_y, radius,
                                   magnification_center=None):
        """Evaluates Gould (2008) eq. 8"""
        shift = radius / sqrt(2.)
        dx = [1., -1., -1., 1.]
        dy = [1., 1., -1., -1.]
        out = []
        for (i, dxval) in enumerate(dx):
            x = source_x + dxval * shift
            y = source_y + dy[i] * shift
            out.append(self.point_source_magnification(
                                              source_x=x, source_y=y))
        if magnification_center is None:
            magnification_center = self.point_source_magnification(
                                    source_x=source_x, source_y=source_y)
        return 0.25 * fsum(out) - magnification_center

    def _rho_check(self, rho):
        """
        Check if rho is float and positive.
        """
        if rho is None:
            raise TypeError(
                'rho must be positive float, but None was provided')
        if not isinstance(rho, float):
            raise TypeError(
                'rho must be positive float, but ' + str(rho) +
                str(type(rho)) + ' was provided')
        if rho < 0:
            raise ValueError(
                'rho must be positive, got: {:}'.format(rho))

    def hexadecapole_magnification(self, source_x, source_y, rho, gamma,
                                   quadrupole=False, all_approximations=False):
        """
        Magnification in hexadecapole approximation of the
        binary-lens/finite-source event - based on `Gould 2008 ApJ
        681, 1593
        <https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract>`_.

        For coordinate system convention see
        :py:func:`point_source_magnification()`

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

            rho: *float*
                Source size relative to Einstein ring radius.

            gamma: *float*
                Linear limb-darkening coefficient in gamma convention.

            quadrupole: *boolean*, optional
                Return quadrupole approximation instead of hexadecapole?
                Default is *False*.

            all_approximations: *boolean*, optional
                Return hexadecapole, quadrupole, and point source
                approximations? Default is *False*.

        Returns :
            magnification: *float* or *sequence* of three *floats*
                Hexadecapole approximation (*float*) by default.
                Quadrupole approximation (*float*) if *quadrupole*
                parameter is *True*. Hexadecapole, quadrupole, and
                point source approximations (*sequence* of three
                *floats*) if *all_approximations* parameter is *True*.
        """
        # In this function, variables named a_* depict magnification.
        if quadrupole and all_approximations:
            raise ValueError('Inconsistent parameters of ' +
                             'BinaryLens.hexadecapole_magnification()')
        self._rho_check(rho)

        a_center = self.point_source_magnification(
            source_x=source_x, source_y=source_y)
        a_rho_half_plus = self._get_magnification_w_plus(
            source_x=source_x, source_y=source_y, radius=0.5*rho,
            magnification_center=a_center)
        a_rho_plus = self._get_magnification_w_plus(
            source_x=source_x, source_y=source_y, radius=rho,
            magnification_center=a_center)

        # This is Gould 2008 eq. 9:
        a_2_rho_square = (16. * a_rho_half_plus - a_rho_plus) / 3.

        # Gould 2008 eq. 6 (part 1/2):
        a_quadrupole = a_center + a_2_rho_square * (1. - 0.2 * gamma)

        # At this point is quadrupole approximation is finished
        if quadrupole:
            return a_quadrupole

        a_rho_times = self._get_magnification_w_times(
            source_x=source_x, source_y=source_y, radius=rho,
            magnification_center=a_center)

        # This is Gould (2008) eq. 9:
        a_4_rho_power4 = 0.5 * (a_rho_plus + a_rho_times) - a_2_rho_square
        # This is Gould (2008) eq. 6 (part 2/2):
        a_add = a_4_rho_power4 * (1. - 11. * gamma / 35.)
        a_hexadecapole = a_quadrupole + a_add

        if all_approximations:
            return (a_hexadecapole, a_quadrupole, a_center)
        else:
            return a_hexadecapole

    def adaptive_contouring_magnification(
            self, source_x, source_y, rho, gamma=None, u_limb_darkening=None,
            accuracy=0.1, ld_accuracy=0.001):
        """
        Binary lens finite source magnification calculated using
        Adaptive Contouring method by `Dominik 2007 MNRAS, 377, 1679
        <https://ui.adsabs.harvard.edu/abs/2007MNRAS.377.1679D/abstract>`_

        See also
        `AdaptiveContouring website by Martin Dominik
        <http://star-www.st-and.ac.uk/~md35/Software.html>`_

        For coordinate system convention see
        :py:func:`point_source_magnification()`

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

            rho: *float*
                Source size relative to Einstein ring radius.

            gamma: *float*, optional
                Linear limb-darkening coefficient in gamma convention.

            u_limb_darkening: *float*
                Linear limb-darkening coefficient in u convention.
                Note that either *gamma* or *u_limb_darkening* can be
                set.  If neither of them is provided then limb
                darkening is ignored.

            accuracy: *float*, optional
                Requested accuracy of the result defined as the sum of
                the area of the squares that determine the contour
                line and the estimated total enclosed area (see sec. 4
                of the paper).  As M. Dominik states: *"this vastly
                overestimates the fractional error, and a suitable
                value should be chosen by testing how its variation
                affects the final results - I recommend starting at
                acc = 0.1."* It significantly affects execution time.

            ld_accuracy: *float*, optional
                Requested limb-darkening accuracy. As M. Dominik
                states: *" Fractional uncertainty for the adaptive
                Simpson integration of the limb-darkening
                profile-related function during application of Green's
                theorem."* It does not add execution time so can be
                set to very small value.

        Returns :
            magnification: *float*
                Magnification.


        """
        if accuracy <= 0.:
            raise ValueError('adaptive_contouring requires accuracy > 0')
        if ld_accuracy <= 0.:
            raise ValueError('adaptive_contouring requires ld_accuracy > 0')
        # Note that this accuracy is not guaranteed.
        self._rho_check(rho)

        if not _adaptive_contouring_wrapped:
            raise ValueError('Adaptive Contouring was not imported properly')

        if gamma is not None and u_limb_darkening is not None:
            raise ValueError(
                'Only one limb darkening parameters can be set' +
                ' in BinaryLens.adaptive_contouring_magnification()')
        elif gamma is not None:
            gamma = float(gamma)
        elif u_limb_darkening is not None:
            gamma = float(mm.Utils.u_to_gamma(u_limb_darkening))
        else:
            gamma = float(0.0)

        s = float(self.separation)
        q = float(self.mass_2 / self.mass_1)
        # AdaptiveContouring uses different coordinates conventions,
        # so we have to transform the coordinates below.
        x = float(-source_x)
        y = float(-source_y)
        rho = float(rho)
        accuracy = float(accuracy)
        ld_accuracy = float(ld_accuracy)

        magnification = _adaptive_contouring_linear(
            s, q, x, y, rho, gamma, accuracy, ld_accuracy)

        return magnification

    def vbbl_magnification(self, source_x, source_y, rho,
                           gamma=None, u_limb_darkening=None,
                           accuracy=0.001):
        """
        Binary lens finite source magnification calculated using VBBL
        library that implements advanced contour integration algorithm
        presented by `Bozza 2010 MNRAS, 408, 2188
        <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract>`_.
        See also `VBBL website by Valerio Bozza
        <http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm>`_.

        For coordinate system convention see
        :py:func:`point_source_magnification()`

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

            rho: *float*
                Source size relative to Einstein ring radius.

            gamma: *float*, optional
                Linear limb-darkening coefficient in gamma convention.

            u_limb_darkening: *float*
                Linear limb-darkening coefficient in u convention.
                Note that either *gamma* or *u_limb_darkening* can be
                set.  If neither of them is provided then limb
                darkening is ignored.

            accuracy: *float*, optional
                Requested accuracy of the result.

        Returns :
            magnification: *float*
                Magnification.

        """
        self._rho_check(rho)
        if accuracy <= 0.:
            raise ValueError(
                "VBBL requires accuracy > 0 e.g. 0.01 or 0.001;" +
                "\n{:} was  provided".format(accuracy))

        if not _vbbl_wrapped:
            raise ValueError('VBBL was not imported properly')

        if gamma is not None and u_limb_darkening is not None:
            raise ValueError('Only one limb darkening parameters can be set' +
                             ' in BinaryLens.vbbl_magnification()')
        elif gamma is not None:
            u_limb_darkening = float(mm.Utils.gamma_to_u(gamma))
        elif u_limb_darkening is not None:
            u_limb_darkening = float(u_limb_darkening)
        else:
            u_limb_darkening = float(0.0)

        s = float(self.separation)
        q = float(self.mass_2 / self.mass_1)
        x = float(source_x)
        y = float(source_y)
        rho = float(rho)
        accuracy = float(accuracy)

        magnification = _vbbl_binary_mag_dark(
            s, q, x, y, rho, u_limb_darkening, accuracy)

        return magnification

class BinaryLensWithShear(BinaryLens):
    """
    The binary lens equation - its solutions, images, parities,
    magnifications, etc.

    The binary lens equation is a 5th order complex polynomial 
    or a 9th order complex polynomial if including external shear.

    Attributes :
        mass_1: *float*
            mass of the primary (left-hand object) as a fraction of
            the total mass.

        mass_2: *float*
            mass of the secondary (right-hand object) as a fraction of the
            total mass.

        separation: *float*
            separation between the two bodies as a fraction of the Einstein
            ring.
        
        convergence_K: *float*
            External mass sheet convergence.
        
        shear_G: *complex*
            External mass sheat shear.

    Note: mass_1 and mass_2 may be defined as a fraction of some other
    mass than the total mass. This is possible but not recommended -
    make sure you know what you're doing before you start using this
    possibility.

    """
    def __init__(self, mass_1=None, mass_2=None, separation=None, convergence_K=0.0, shear_G=complex(0,0)):
        BinaryLens.__init__(self, mass_1, mass_2, separation)
        self.convergence_K = convergence_K
        self.shear_G = shear_G

    def _get_polynomial_WM95(self, source_x, source_y):
        """
        calculate coefficients of the polynomial in geometric center frame
        """
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff

        pos_z1 = self._position_z1

        z1_pow2 = pos_z1 * pos_z1
        z1_pow3 = z1_pow2 * pos_z1
        z1_pow4 = z1_pow2 * z1_pow2
        z1_pow5 = z1_pow2 * z1_pow2 * pos_z1
        z1_pow6 = z1_pow2 * z1_pow2 * z1_pow2

        # #Convergence added here
        convergence_K = self.convergence_K
        K_pow2 = convergence_K * convergence_K
        K_pow3 = convergence_K * K_pow2

        zeta = self._zeta
        zeta_conj = zeta.conjugate()
        zeta_conj_pow2 = zeta_conj * zeta_conj

        # # Calculate the coefficients of the 5th order complex polynomial now with external convergence
        coeff_5 = mm.Utils.complex_fsum([-z1_pow2*K_pow3 , 3*z1_pow2*K_pow2 , - 3*z1_pow2*convergence_K , z1_pow2 , convergence_K*zeta_conj_pow2 , - zeta_conj_pow2])
        coeff_4 = mm.Utils.complex_fsum([2*total_m*convergence_K*zeta_conj , - 2*total_m*zeta_conj , - zeta*z1_pow2*K_pow2 , 2*zeta*z1_pow2*convergence_K , - zeta*z1_pow2 , zeta*zeta_conj_pow2 , 2*pos_z1*K_pow2*m_diff , - 4*pos_z1*convergence_K*m_diff , 2*pos_z1*m_diff])
        coeff_3 = mm.Utils.complex_fsum([4*total_m*zeta*zeta_conj , 2*z1_pow4*K_pow3 , - 6*z1_pow4*K_pow2 , 6*z1_pow4*convergence_K , - 2*z1_pow4 , - 2*z1_pow2*convergence_K*zeta_conj_pow2 , 2*z1_pow2*zeta_conj_pow2 , 4*pos_z1*convergence_K*m_diff*zeta_conj , - 4*pos_z1*m_diff*zeta_conj])
        coeff_2 = mm.Utils.complex_fsum([4*total_m_pow2*zeta , 4*total_m*pos_z1*convergence_K*m_diff , - 4*total_m*pos_z1*m_diff , 2*zeta*z1_pow4*K_pow2 , - 4*zeta*z1_pow4*convergence_K , 2*zeta*z1_pow4 , - 2*zeta*z1_pow2*zeta_conj_pow2 , 4*zeta*pos_z1*m_diff*zeta_conj , - 4*z1_pow3*K_pow2*m_diff , 8*z1_pow3*convergence_K*m_diff , - 4*z1_pow3*m_diff])
        coeff_1 = mm.Utils.complex_fsum([4*total_m_pow2*z1_pow2*convergence_K , - 4*total_m_pow2*z1_pow2 , - 4*total_m*zeta*z1_pow2*zeta_conj , 8*total_m*zeta*pos_z1*m_diff , - z1_pow6*K_pow3 , 3*z1_pow6*K_pow2 , - 3*z1_pow6*convergence_K , z1_pow6 , z1_pow4*convergence_K*zeta_conj_pow2 , - z1_pow4*zeta_conj_pow2 , - 4*z1_pow3*convergence_K*m_diff*zeta_conj , 4*z1_pow3*m_diff*zeta_conj , 4*z1_pow2*convergence_K*m_diff_pow2 , - 4*z1_pow2*m_diff_pow2])
        coeff_0 = mm.Utils.complex_fsum([-2*total_m*z1_pow4*convergence_K*zeta_conj , 2*total_m*z1_pow4*zeta_conj , 4*total_m*z1_pow3*convergence_K*m_diff , - 4*total_m*z1_pow3*m_diff , - zeta*z1_pow6*K_pow2 , 2*zeta*z1_pow6*convergence_K , - zeta*z1_pow6 , zeta*z1_pow4*zeta_conj_pow2 , - 4*zeta*z1_pow3*m_diff*zeta_conj , 4*zeta*z1_pow2*m_diff_pow2 , 2*z1_pow5*K_pow2*m_diff , - 4*z1_pow5*convergence_K*m_diff , 2*z1_pow5*m_diff])

        # Return the coefficients of the polynomial
        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5]
        return np.array(coeffs_list).reshape(6)

    def _get_polynomial_planet_frame(self, source_x, source_y):
        """calculate coefficients of the polynomial in planet frame"""
        # Calculate constants
        self._calculate_variables(source_x=source_x, source_y=source_y)
        total_m = self._total_mass
        total_m_pow2 = total_m * total_m
        total_m_pow3 = total_m * total_m_pow2

        m_diff = self._mass_difference
        m_diff_pow2 = m_diff * m_diff
        m_diff_pow3 = m_diff * m_diff_pow2

        zeta = self._zeta
        zeta_conj = zeta.conjugate()
        zeta_conj_pow2 = zeta_conj * zeta_conj
        zeta_conj_pow3 = zeta_conj * zeta_conj_pow2

        # External Convergence added here
        convergence_K = self.convergence_K
        K_pow2 = convergence_K * convergence_K
        K_pow3 = convergence_K * K_pow2
        K_pow4 = K_pow2 * K_pow2

        #External Shear added here
        shear_G = self.shear_G
        Gc = np.conjugate(shear_G)
        G_pow2 = shear_G * shear_G
        G_pow3 = shear_G * G_pow2
        G_pow4 = G_pow2 * G_pow2

        Gc_pow2 = Gc * Gc
        Gc_pow3 = Gc * Gc_pow2
        Gc_pow4 = Gc_pow2 * Gc_pow2
    
        c_sum = mm.Utils.complex_fsum

        z1 = self._position_z1
        z1_pow2 = z1 * z1
        z1_pow3 = z1_pow2 * z1
        z1_pow4 = z1_pow2 * z1_pow2

        coeff_9 = c_sum([-shear_G*Gc_pow3 , Gc_pow2*K_pow2 , - 2*Gc_pow2*convergence_K , Gc_pow2])
        coeff_8 = c_sum([zeta*Gc_pow2*convergence_K , - zeta*Gc_pow2 , 3*z1*shear_G*Gc_pow3 , - z1*shear_G*Gc_pow2*convergence_K , z1*shear_G*Gc_pow2 , - 3*z1*Gc_pow2*K_pow2 , 6*z1*Gc_pow2*convergence_K , - 3*z1*Gc_pow2 , z1*Gc*K_pow3 , - 3*z1*Gc*K_pow2 , 3*z1*Gc*convergence_K , - z1*Gc , - 3*shear_G*Gc_pow2*zeta_conj , 2*Gc*K_pow2*zeta_conj , - 4*Gc*convergence_K*zeta_conj , 2*Gc*zeta_conj])
        coeff_7 = c_sum([-6*total_m*shear_G*Gc_pow2 , 2*total_m*Gc*K_pow2 , - 4*total_m*Gc*convergence_K , 2*total_m*Gc , - 3*zeta*z1*Gc_pow2*convergence_K , 3*zeta*z1*Gc_pow2 , zeta*z1*Gc*K_pow2 , - 2*zeta*z1*Gc*convergence_K , zeta*z1*Gc , 2*zeta*Gc*convergence_K*zeta_conj , - 2*zeta*Gc*zeta_conj , - 3*z1_pow2*shear_G*Gc_pow3 , 3*z1_pow2*shear_G*Gc_pow2*convergence_K , - 3*z1_pow2*shear_G*Gc_pow2 , 3*z1_pow2*Gc_pow2*K_pow2 , - 6*z1_pow2*Gc_pow2*convergence_K , 3*z1_pow2*Gc_pow2 , - 3*z1_pow2*Gc*K_pow3 , 9*z1_pow2*Gc*K_pow2 , - 9*z1_pow2*Gc*convergence_K , 3*z1_pow2*Gc , 9*z1*shear_G*Gc_pow2*zeta_conj , - 2*z1*shear_G*Gc*convergence_K*zeta_conj , 2*z1*shear_G*Gc*zeta_conj , - 6*z1*Gc*K_pow2*zeta_conj , 12*z1*Gc*convergence_K*zeta_conj , - 6*z1*Gc*zeta_conj , z1*K_pow3*zeta_conj , - 3*z1*K_pow2*zeta_conj , 3*z1*convergence_K*zeta_conj , - z1*zeta_conj , - 3*shear_G*Gc*zeta_conj_pow2 , K_pow2*zeta_conj_pow2 , - 2*convergence_K*zeta_conj_pow2 , zeta_conj_pow2])
        coeff_6 = c_sum([4*total_m*zeta*Gc*convergence_K , - 4*total_m*zeta*Gc , 15*total_m*z1*shear_G*Gc_pow2 , - 4*total_m*z1*shear_G*Gc*convergence_K , 4*total_m*z1*shear_G*Gc , - 4*total_m*z1*Gc*K_pow2 , 8*total_m*z1*Gc*convergence_K , - 4*total_m*z1*Gc , total_m*z1*K_pow3 , - 3*total_m*z1*K_pow2 , 3*total_m*z1*convergence_K , - total_m*z1 , - 12*total_m*shear_G*Gc*zeta_conj , 2*total_m*K_pow2*zeta_conj , - 4*total_m*convergence_K*zeta_conj , 2*total_m*zeta_conj , 3*zeta*z1_pow2*Gc_pow2*convergence_K , - 3*zeta*z1_pow2*Gc_pow2 , - 3*zeta*z1_pow2*Gc*K_pow2 , 6*zeta*z1_pow2*Gc*convergence_K , - 3*zeta*z1_pow2*Gc , - 6*zeta*z1*Gc*convergence_K*zeta_conj , 6*zeta*z1*Gc*zeta_conj , zeta*z1*K_pow2*zeta_conj , - 2*zeta*z1*convergence_K*zeta_conj , zeta*z1*zeta_conj , zeta*convergence_K*zeta_conj_pow2 , - zeta*zeta_conj_pow2 , z1_pow3*shear_G*Gc_pow3 , - 3*z1_pow3*shear_G*Gc_pow2*convergence_K , 3*z1_pow3*shear_G*Gc_pow2 , - z1_pow3*Gc_pow2*K_pow2 , 2*z1_pow3*Gc_pow2*convergence_K , - z1_pow3*Gc_pow2 , 3*z1_pow3*Gc*K_pow3 , - 9*z1_pow3*Gc*K_pow2 , 9*z1_pow3*Gc*convergence_K , - 3*z1_pow3*Gc , - 9*z1_pow2*shear_G*Gc_pow2*zeta_conj , 6*z1_pow2*shear_G*Gc*convergence_K*zeta_conj , - 6*z1_pow2*shear_G*Gc*zeta_conj , 6*z1_pow2*Gc*K_pow2*zeta_conj , - 12*z1_pow2*Gc*convergence_K*zeta_conj , 6*z1_pow2*Gc*zeta_conj , - 3*z1_pow2*K_pow3*zeta_conj , 9*z1_pow2*K_pow2*zeta_conj , - 9*z1_pow2*convergence_K*zeta_conj , 3*z1_pow2*zeta_conj , 3*z1*shear_G*Gc_pow2*m_diff , 9*z1*shear_G*Gc*zeta_conj_pow2 , - z1*shear_G*convergence_K*zeta_conj_pow2 , z1*shear_G*zeta_conj_pow2 , - 2*z1*Gc*K_pow2*m_diff , 4*z1*Gc*convergence_K*m_diff , - 2*z1*Gc*m_diff , - z1*K_pow3*m_diff , 3*z1*K_pow2*m_diff , - 3*z1*K_pow2*zeta_conj_pow2 , - 3*z1*convergence_K*m_diff , 6*z1*convergence_K*zeta_conj_pow2 , z1*m_diff , - 3*z1*zeta_conj_pow2 , - shear_G*zeta_conj_pow3])
        coeff_5 = c_sum([-12*total_m_pow2*shear_G*Gc , - 10*total_m*zeta*z1*Gc*convergence_K , 10*total_m*zeta*z1*Gc , 2*total_m*zeta*z1*K_pow2 , - 4*total_m*zeta*z1*convergence_K , 2*total_m*zeta*z1 , 4*total_m*zeta*convergence_K*zeta_conj , - 4*total_m*zeta*zeta_conj , - 12*total_m*z1_pow2*shear_G*Gc_pow2 , 10*total_m*z1_pow2*shear_G*Gc*convergence_K , - 10*total_m*z1_pow2*shear_G*Gc , 2*total_m*z1_pow2*Gc*K_pow2 , - 4*total_m*z1_pow2*Gc*convergence_K , 2*total_m*z1_pow2*Gc , - 2*total_m*z1_pow2*K_pow3 , 6*total_m*z1_pow2*K_pow2 , - 6*total_m*z1_pow2*convergence_K , 2*total_m*z1_pow2 , 30*total_m*z1*shear_G*Gc*zeta_conj , - 4*total_m*z1*shear_G*convergence_K*zeta_conj , 4*total_m*z1*shear_G*zeta_conj , - 4*total_m*z1*K_pow2*zeta_conj , 8*total_m*z1*convergence_K*zeta_conj , - 4*total_m*z1*zeta_conj , - 6*total_m*shear_G*zeta_conj_pow2 , - zeta*z1_pow3*Gc_pow2*convergence_K , zeta*z1_pow3*Gc_pow2 , 3*zeta*z1_pow3*Gc*K_pow2 , - 6*zeta*z1_pow3*Gc*convergence_K , 3*zeta*z1_pow3*Gc , 6*zeta*z1_pow2*Gc*convergence_K*zeta_conj , - 6*zeta*z1_pow2*Gc*zeta_conj , - 3*zeta*z1_pow2*K_pow2*zeta_conj , 6*zeta*z1_pow2*convergence_K*zeta_conj , - 3*zeta*z1_pow2*zeta_conj , - 2*zeta*z1*Gc*convergence_K*m_diff , 2*zeta*z1*Gc*m_diff , - 3*zeta*z1*convergence_K*zeta_conj_pow2 , 3*zeta*z1*zeta_conj_pow2 , z1_pow4*shear_G*Gc_pow2*convergence_K , - z1_pow4*shear_G*Gc_pow2 , - z1_pow4*Gc*K_pow3 , 3*z1_pow4*Gc*K_pow2 , - 3*z1_pow4*Gc*convergence_K , z1_pow4*Gc , 3*z1_pow3*shear_G*Gc_pow2*zeta_conj , - 6*z1_pow3*shear_G*Gc*convergence_K*zeta_conj , 6*z1_pow3*shear_G*Gc*zeta_conj , - 2*z1_pow3*Gc*K_pow2*zeta_conj , 4*z1_pow3*Gc*convergence_K*zeta_conj , - 2*z1_pow3*Gc*zeta_conj , 3*z1_pow3*K_pow3*zeta_conj , - 9*z1_pow3*K_pow2*zeta_conj , 9*z1_pow3*convergence_K*zeta_conj , - 3*z1_pow3*zeta_conj , - 6*z1_pow2*shear_G*Gc_pow2*m_diff , 2*z1_pow2*shear_G*Gc*convergence_K*m_diff , - 2*z1_pow2*shear_G*Gc*m_diff , - 9*z1_pow2*shear_G*Gc*zeta_conj_pow2 , 3*z1_pow2*shear_G*convergence_K*zeta_conj_pow2 , - 3*z1_pow2*shear_G*zeta_conj_pow2 , 4*z1_pow2*Gc*K_pow2*m_diff , - 8*z1_pow2*Gc*convergence_K*m_diff , 4*z1_pow2*Gc*m_diff , 2*z1_pow2*K_pow3*m_diff , - 6*z1_pow2*K_pow2*m_diff , 3*z1_pow2*K_pow2*zeta_conj_pow2 , 6*z1_pow2*convergence_K*m_diff , - 6*z1_pow2*convergence_K*zeta_conj_pow2 , - 2*z1_pow2*m_diff , 3*z1_pow2*zeta_conj_pow2 , 6*z1*shear_G*Gc*m_diff*zeta_conj , 3*z1*shear_G*zeta_conj_pow3 , - 2*z1*K_pow2*m_diff*zeta_conj , 4*z1*convergence_K*m_diff*zeta_conj , - 2*z1*m_diff*zeta_conj])
        coeff_4 = c_sum([4*total_m_pow2*zeta*convergence_K , - 4*total_m_pow2*zeta , 24*total_m_pow2*z1*shear_G*Gc , - 4*total_m_pow2*z1*shear_G*convergence_K , 4*total_m_pow2*z1*shear_G , 2*total_m_pow2*z1*K_pow2 , - 4*total_m_pow2*z1*convergence_K , 2*total_m_pow2*z1 , - 12*total_m_pow2*shear_G*zeta_conj , 8*total_m*zeta*z1_pow2*Gc*convergence_K , - 8*total_m*zeta*z1_pow2*Gc , - 5*total_m*zeta*z1_pow2*K_pow2 , 10*total_m*zeta*z1_pow2*convergence_K , - 5*total_m*zeta*z1_pow2 , - 10*total_m*zeta*z1*convergence_K*zeta_conj , 10*total_m*zeta*z1*zeta_conj , 3*total_m*z1_pow3*shear_G*Gc_pow2 , - 8*total_m*z1_pow3*shear_G*Gc*convergence_K , 8*total_m*z1_pow3*shear_G*Gc , total_m*z1_pow3*K_pow3 , - 3*total_m*z1_pow3*K_pow2 , 3*total_m*z1_pow3*convergence_K , - total_m*z1_pow3 , - 24*total_m*z1_pow2*shear_G*Gc*zeta_conj , 10*total_m*z1_pow2*shear_G*convergence_K*zeta_conj , - 10*total_m*z1_pow2*shear_G*zeta_conj , 2*total_m*z1_pow2*K_pow2*zeta_conj , - 4*total_m*z1_pow2*convergence_K*zeta_conj , 2*total_m*z1_pow2*zeta_conj , 12*total_m*z1*shear_G*Gc*m_diff , 15*total_m*z1*shear_G*zeta_conj_pow2 , - 2*total_m*z1*K_pow2*m_diff , 4*total_m*z1*convergence_K*m_diff , - 2*total_m*z1*m_diff , - zeta*z1_pow4*Gc*K_pow2 , 2*zeta*z1_pow4*Gc*convergence_K , - zeta*z1_pow4*Gc , - 2*zeta*z1_pow3*Gc*convergence_K*zeta_conj , 2*zeta*z1_pow3*Gc*zeta_conj , 3*zeta*z1_pow3*K_pow2*zeta_conj , - 6*zeta*z1_pow3*convergence_K*zeta_conj , 3*zeta*z1_pow3*zeta_conj , 4*zeta*z1_pow2*Gc*convergence_K*m_diff , - 4*zeta*z1_pow2*Gc*m_diff , - zeta*z1_pow2*K_pow2*m_diff , 2*zeta*z1_pow2*convergence_K*m_diff , 3*zeta*z1_pow2*convergence_K*zeta_conj_pow2 , - zeta*z1_pow2*m_diff , - 3*zeta*z1_pow2*zeta_conj_pow2 , - 2*zeta*z1*convergence_K*m_diff*zeta_conj , 2*zeta*z1*m_diff*zeta_conj , 2*z1_pow4*shear_G*Gc*convergence_K*zeta_conj , - 2*z1_pow4*shear_G*Gc*zeta_conj , - z1_pow4*K_pow3*zeta_conj , 3*z1_pow4*K_pow2*zeta_conj , - 3*z1_pow4*convergence_K*zeta_conj , z1_pow4*zeta_conj , 3*z1_pow3*shear_G*Gc_pow2*m_diff , - 4*z1_pow3*shear_G*Gc*convergence_K*m_diff , 4*z1_pow3*shear_G*Gc*m_diff , 3*z1_pow3*shear_G*Gc*zeta_conj_pow2 , - 3*z1_pow3*shear_G*convergence_K*zeta_conj_pow2 , 3*z1_pow3*shear_G*zeta_conj_pow2 , - 2*z1_pow3*Gc*K_pow2*m_diff , 4*z1_pow3*Gc*convergence_K*m_diff , - 2*z1_pow3*Gc*m_diff , - z1_pow3*K_pow3*m_diff , 3*z1_pow3*K_pow2*m_diff , - z1_pow3*K_pow2*zeta_conj_pow2 , - 3*z1_pow3*convergence_K*m_diff , 2*z1_pow3*convergence_K*zeta_conj_pow2 , z1_pow3*m_diff , - z1_pow3*zeta_conj_pow2 , - 12*z1_pow2*shear_G*Gc*m_diff*zeta_conj , 2*z1_pow2*shear_G*convergence_K*m_diff*zeta_conj , - 2*z1_pow2*shear_G*m_diff*zeta_conj , - 3*z1_pow2*shear_G*zeta_conj_pow3 , 4*z1_pow2*K_pow2*m_diff*zeta_conj , - 8*z1_pow2*convergence_K*m_diff*zeta_conj , 4*z1_pow2*m_diff*zeta_conj , 3*z1*shear_G*m_diff*zeta_conj_pow2])
        coeff_3 = c_sum([-8*total_m_pow3*shear_G , - 8*total_m_pow2*zeta*z1*convergence_K , 8*total_m_pow2*zeta*z1 , - 15*total_m_pow2*z1_pow2*shear_G*Gc , 8*total_m_pow2*z1_pow2*shear_G*convergence_K , - 8*total_m_pow2*z1_pow2*shear_G , - 3*total_m_pow2*z1_pow2*K_pow2 , 6*total_m_pow2*z1_pow2*convergence_K , - 3*total_m_pow2*z1_pow2 , 24*total_m_pow2*z1*shear_G*zeta_conj , - 2*total_m*zeta*z1_pow3*Gc*convergence_K , 2*total_m*zeta*z1_pow3*Gc , 4*total_m*zeta*z1_pow3*K_pow2 , - 8*total_m*zeta*z1_pow3*convergence_K , 4*total_m*zeta*z1_pow3 , 8*total_m*zeta*z1_pow2*convergence_K*zeta_conj , - 8*total_m*zeta*z1_pow2*zeta_conj , - 4*total_m*zeta*z1*convergence_K*m_diff , 4*total_m*zeta*z1*m_diff , 2*total_m*z1_pow4*shear_G*Gc*convergence_K , - 2*total_m*z1_pow4*shear_G*Gc , 6*total_m*z1_pow3*shear_G*Gc*zeta_conj , - 8*total_m*z1_pow3*shear_G*convergence_K*zeta_conj , 8*total_m*z1_pow3*shear_G*zeta_conj , - 18*total_m*z1_pow2*shear_G*Gc*m_diff , 4*total_m*z1_pow2*shear_G*convergence_K*m_diff , - 4*total_m*z1_pow2*shear_G*m_diff , - 12*total_m*z1_pow2*shear_G*zeta_conj_pow2 , 2*total_m*z1_pow2*K_pow2*m_diff , - 4*total_m*z1_pow2*convergence_K*m_diff , 2*total_m*z1_pow2*m_diff , 12*total_m*z1*shear_G*m_diff*zeta_conj , - zeta*z1_pow4*K_pow2*zeta_conj , 2*zeta*z1_pow4*convergence_K*zeta_conj , - zeta*z1_pow4*zeta_conj , - 2*zeta*z1_pow3*Gc*convergence_K*m_diff , 2*zeta*z1_pow3*Gc*m_diff , 2*zeta*z1_pow3*K_pow2*m_diff , - 4*zeta*z1_pow3*convergence_K*m_diff , - zeta*z1_pow3*convergence_K*zeta_conj_pow2 , 2*zeta*z1_pow3*m_diff , zeta*z1_pow3*zeta_conj_pow2 , 4*zeta*z1_pow2*convergence_K*m_diff*zeta_conj , - 4*zeta*z1_pow2*m_diff*zeta_conj , 2*z1_pow4*shear_G*Gc*convergence_K*m_diff , - 2*z1_pow4*shear_G*Gc*m_diff , z1_pow4*shear_G*convergence_K*zeta_conj_pow2 , - z1_pow4*shear_G*zeta_conj_pow2 , 6*z1_pow3*shear_G*Gc*m_diff*zeta_conj , - 4*z1_pow3*shear_G*convergence_K*m_diff*zeta_conj , 4*z1_pow3*shear_G*m_diff*zeta_conj , z1_pow3*shear_G*zeta_conj_pow3 , - 2*z1_pow3*K_pow2*m_diff*zeta_conj , 4*z1_pow3*convergence_K*m_diff*zeta_conj , - 2*z1_pow3*m_diff*zeta_conj , - 3*z1_pow2*shear_G*Gc*m_diff_pow2 , - 6*z1_pow2*shear_G*m_diff*zeta_conj_pow2 , z1_pow2*K_pow2*m_diff_pow2 , - 2*z1_pow2*convergence_K*m_diff_pow2 , z1_pow2*m_diff_pow2])
        coeff_2 = c_sum([12*total_m_pow3*z1*shear_G , 5*total_m_pow2*zeta*z1_pow2*convergence_K , - 5*total_m_pow2*zeta*z1_pow2 , 3*total_m_pow2*z1_pow3*shear_G*Gc , - 5*total_m_pow2*z1_pow3*shear_G*convergence_K , 5*total_m_pow2*z1_pow3*shear_G , total_m_pow2*z1_pow3*K_pow2 , - 2*total_m_pow2*z1_pow3*convergence_K , total_m_pow2*z1_pow3 , - 15*total_m_pow2*z1_pow2*shear_G*zeta_conj , 12*total_m_pow2*z1*shear_G*m_diff , - total_m*zeta*z1_pow4*K_pow2 , 2*total_m*zeta*z1_pow4*convergence_K , - total_m*zeta*z1_pow4 , - 2*total_m*zeta*z1_pow3*convergence_K*zeta_conj , 2*total_m*zeta*z1_pow3*zeta_conj , 6*total_m*zeta*z1_pow2*convergence_K*m_diff , - 6*total_m*zeta*z1_pow2*m_diff , 2*total_m*z1_pow4*shear_G*convergence_K*zeta_conj , - 2*total_m*z1_pow4*shear_G*zeta_conj , 6*total_m*z1_pow3*shear_G*Gc*m_diff , - 6*total_m*z1_pow3*shear_G*convergence_K*m_diff , 6*total_m*z1_pow3*shear_G*m_diff , 3*total_m*z1_pow3*shear_G*zeta_conj_pow2 , - 18*total_m*z1_pow2*shear_G*m_diff*zeta_conj , - zeta*z1_pow4*K_pow2*m_diff , 2*zeta*z1_pow4*convergence_K*m_diff , - zeta*z1_pow4*m_diff , - 2*zeta*z1_pow3*convergence_K*m_diff*zeta_conj , 2*zeta*z1_pow3*m_diff*zeta_conj , zeta*z1_pow2*convergence_K*m_diff_pow2 , - zeta*z1_pow2*m_diff_pow2 , 2*z1_pow4*shear_G*convergence_K*m_diff*zeta_conj , - 2*z1_pow4*shear_G*m_diff*zeta_conj , 3*z1_pow3*shear_G*Gc*m_diff_pow2 , - z1_pow3*shear_G*convergence_K*m_diff_pow2 , z1_pow3*shear_G*m_diff_pow2 , 3*z1_pow3*shear_G*m_diff*zeta_conj_pow2 , - z1_pow3*K_pow2*m_diff_pow2 , 2*z1_pow3*convergence_K*m_diff_pow2 , - z1_pow3*m_diff_pow2 , - 3*z1_pow2*shear_G*m_diff_pow2*zeta_conj])
        coeff_1 = c_sum([-6*total_m_pow3*z1_pow2*shear_G , - total_m_pow2*zeta*z1_pow3*convergence_K , total_m_pow2*zeta*z1_pow3 , total_m_pow2*z1_pow4*shear_G*convergence_K , - total_m_pow2*z1_pow4*shear_G , 3*total_m_pow2*z1_pow3*shear_G*zeta_conj , - 12*total_m_pow2*z1_pow2*shear_G*m_diff , - 2*total_m*zeta*z1_pow3*convergence_K*m_diff , 2*total_m*zeta*z1_pow3*m_diff , 2*total_m*z1_pow4*shear_G*convergence_K*m_diff , - 2*total_m*z1_pow4*shear_G*m_diff , 6*total_m*z1_pow3*shear_G*m_diff*zeta_conj , - 6*total_m*z1_pow2*shear_G*m_diff_pow2 , - zeta*z1_pow3*convergence_K*m_diff_pow2 , zeta*z1_pow3*m_diff_pow2 , z1_pow4*shear_G*convergence_K*m_diff_pow2 , - z1_pow4*shear_G*m_diff_pow2 , 3*z1_pow3*shear_G*m_diff_pow2*zeta_conj])
        coeff_0 = c_sum([total_m_pow3*z1_pow3*shear_G , 3*total_m_pow2*z1_pow3*shear_G*m_diff , 3*total_m*z1_pow3*shear_G*m_diff_pow2 , z1_pow3*shear_G*m_diff_pow3])

        coeffs_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_6, coeff_7, coeff_8, coeff_9]
        return np.array(coeffs_list).reshape(10)

    def _get_polynomial_roots(self, source_x, source_y):
        """roots of the polynomial"""
        polynomial_input = [self.mass_1, self.mass_2, self.separation, self.convergence_K, self.shear_G,
                            source_x, source_y]

        if polynomial_input == self._last_polynomial_input:
            return self._polynomial_roots

        polynomial = self._get_polynomial(
            source_x=source_x, source_y=source_y)

        np_polyroots = np.polynomial.polynomial.polyroots
        if self._solver == 'numpy':
            self._polynomial_roots = np_polyroots(polynomial)
        elif self._solver == 'Skowron_and_Gould_12':
            args = polynomial.real.tolist() + polynomial.imag.tolist()
            try:
                out = _vbbl_SG12_9(*args)
            except ValueError as err:
                err2 = "\n\nSwitching from Skowron & Gould 2012 to numpy"
                warnings.warn(str(err) + err2, UserWarning)
                self._solver = 'numpy'
                self._polynomial_roots = np_polyroots(polynomial)
            else:
                roots = [out[i] + out[i+9] * 1.j for i in range(9) if ((abs(out[i]) > 1e-10 or abs(out[i+9]) > 1e-10) 
                        and (abs(out[i] - self._position_z1) > 1e-10 or abs(out[i+9]) > 1e-10))]
                self._polynomial_roots = np.array(roots)
        else:
            raise ValueError('Unknown solver: {:}'.format(self._solver))
        self._last_polynomial_input = polynomial_input

        return self._polynomial_roots

    def _polynomial_roots_ok(
            self, source_x, source_y, return_distances=False):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots(
            source_x=source_x, source_y=source_y)

        roots_conj = np.conjugate(roots)
        component2 = self.mass_1 / (roots_conj - self._position_z1)
        component3 = self.mass_2 / (roots_conj - self._position_z2)
        solutions = (self._zeta + self.shear_G * roots_conj + component2 + component3) / (1 - self.convergence_K)

        # This backs-up the lens equation.


        out = []
        distances = []
        for (i, root) in enumerate(roots):
            distances_from_root = abs((solutions-root)**2)
            min_distance_arg = np.argmin(distances_from_root)

            if i == min_distance_arg:
                out.append(root)
                distances.append(distances_from_root[min_distance_arg])
            # The values in distances[] are a diagnostic on how good the
            # numerical accuracy is.

        # If the lens equation is solved correctly, there should be
        # either 3 or 5 solutions (corresponding to 3 or 5 images)
        if len(out) not in [1,2,3,4,5,6,7,8,9]:
            msg = ("Wrong number of solutions to the lens equation of binary" +
                   " lens.\nGot {:} and expected 3 or 5.\nThe parameters " +
                   "(m1, m2, s, source_x, source_y, solver) are:\n" +
                   "{:} {:} {:} {:} {:}  {:}\n\n" +
                   "Consider using 'point_source_point_lens' method for " +
                   "epochs when the source is very far from the lens. Note " +
                   "that it's different from 'point_source' method.")
            txt = msg.format(
                len(out), repr(self.mass_1), repr(self.mass_2),
                repr(self.separation), repr(source_x), repr(source_y),
                self._solver)

            if self._solver != "Skowron_and_Gould_12":
                txt += (
                    "\n\nYou should switch to using Skowron_and_Gould_12" +
                    " polynomial root solver. It is much more accurate than " +
                    "numpy.polynomial.polynomial.polyroots(). " +
                    "Skowron_and_Gould_12 method is selected in automated " +
                    "way if VBBL is imported properly.")
            distance = sqrt(source_x**2 + source_y**2)
            if (self.mass_2 > 1.e-6 * self.mass_1 and
                    (distance < 15. or distance < 2. * self.separation)):
                txt += ("\n\nThis is surprising error - please contact code " +
                        "authors and provide the above error message.")
            elif distance > 200.:
                txt += ("\n\nYou try to calculate magnification at huge " +
                        "distance from the source and this is causing an " +
                        "error.")
            txt += "\nMulensModel version: {:}".format(mm.__version__)

            raise ValueError(txt)
        if return_distances:
            return (np.array(out), np.array(distances))
        else:
            return np.array(out)

    def _jacobian_determinant_ok(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._polynomial_roots_ok(
                                   source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        add_1 = self.mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self.mass_2 / (self._position_z2 - roots_ok_bar)**2
        derivative = add_1 + add_2 - self.shear_G

        return (1.- self.convergence_K)**2 - derivative * np.conjugate(derivative)

    def vbbl_magnification(self, source_x, source_y, rho,
                           gamma=None, u_limb_darkening=None,
                           accuracy=0.001):
        """
        Binary lens finite source magnification calculated using VBBL
        library that implements advanced contour integration algorithm
        presented by `Bozza 2010 MNRAS, 408, 2188
        <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract>`_.
        See also `VBBL website by Valerio Bozza
        <http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm>`_.

        For coordinate system convention see
        :py:func:`point_source_magnification()`

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

            rho: *float*
                Source size relative to Einstein ring radius.

            gamma: *float*, optional
                Linear limb-darkening coefficient in gamma convention.

            u_limb_darkening: *float*
                Linear limb-darkening coefficient in u convention.
                Note that either *gamma* or *u_limb_darkening* can be
                set.  If neither of them is provided then limb
                darkening is ignored.

            accuracy: *float*, optional
                Requested accuracy of the result.

        Returns :
            magnification: *float*
                Magnification.

        """
        self._rho_check(rho)
        if accuracy <= 0.:
            raise ValueError(
                "VBBL requires accuracy > 0 e.g. 0.01 or 0.001;" +
                "\n{:} was  provided".format(accuracy))
        
        if not _vbbl_wrapped:
            raise ValueError('VBBL was not imported properly')

        if gamma is not None and u_limb_darkening is not None:
            raise ValueError('Only one limb darkening parameters can be set' +
                             ' in BinaryLens.vbbl_magnification()')
        elif gamma is not None:
            u_limb_darkening = float(mm.Utils.gamma_to_u(gamma))
        elif u_limb_darkening is not None:
            u_limb_darkening = float(u_limb_darkening)
        else:
            u_limb_darkening = float(0.0)

        s = float(self.separation)
        q = float(self.mass_2 / self.mass_1)
        x = float(source_x)
        y = float(source_y)
        rho = float(rho)
        accuracy = float(accuracy)

        magnification = _vbbl_binary_mag_0(
            s, q, x, y, self.convergence_K, self.shear_G.real, self.shear_G.imag)

        return magnification
