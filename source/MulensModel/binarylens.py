import warnings
import numpy as np
from math import fsum, sqrt

from MulensModel.binarylensimports import (
    _vbbl_wrapped, _adaptive_contouring_wrapped, _vbbl_binary_mag_dark,
    _vbbl_SG12_5, _adaptive_contouring_linear, _solver)
import MulensModel as mm


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

    def _verify_polynomial_roots(
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

    def _get_jacobian_determinant(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._verify_polynomial_roots(
            source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        add_1 = self.mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self.mass_2 / (self._position_z2 - roots_ok_bar)**2
        derivative = add_1 + add_2

        return 1. - derivative * np.conjugate(derivative)

    def _get_signed_magnification(self, source_x, source_y):
        """signed magnification for each image separately"""
        return 1. / self._get_jacobian_determinant(
            source_x=source_x, source_y=source_y)

    def _get_point_source(self, source_x, source_y):
        """calculate point source magnification"""
        signed_magnification = self._get_signed_magnification(
            source_x=source_x, source_y=source_y)
        return fsum(abs(signed_magnification))

    def _get_point_source_Witt_Mao_95(self, source_x, source_y):
        """calculate point source magnification"""
        return self._get_point_source(source_x=source_x, source_y=source_y)

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
        repeat = False
        try:
            out = self._point_source_magnification_VBBL(source_x, source_y)
        except Exception:
            repeat = True

        if repeat or out < 1.:
            out = self._point_source_magnification(source_x, source_y)
        return out

    def _point_source_magnification_VBBL(self, source_x, source_y):
        """
        Calculate point source magnification using VBBL fully
        """
        args = [self.separation, self.mass_2/self.mass_1, source_x, source_y]
        return _vbbl_binary_mag(*[float(arg) for arg in args])

    def _point_source_magnification(self, source_x, source_y):
        """
        Calculate point source magnification using VBBL for solving
        the polynomial and MM code for rest.
        """
        if self._use_planet_frame:
            x_shift = -self.mass_1 / (self.mass_1 + self.mass_2)
        else:
            x_shift = self.mass_2 / (self.mass_1 + self.mass_2) - 0.5
        x_shift *= self.separation
        # We need to add this because in order to shift to correct frame.
        return self._get_point_source_Witt_Mao_95(
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

    def _check_rho(self, rho):
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
        self._check_rho(rho)

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
        self._check_rho(rho)

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
        self._check_rho(rho)
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
