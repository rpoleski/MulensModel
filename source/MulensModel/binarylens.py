import warnings
import numpy as np
from math import fsum, sqrt

from MulensModel.binarylensimports import (
    _vbbl_wrapped, _adaptive_contouring_wrapped,
    _vbbl_binary_mag_dark, _vbbl_binary_mag_finite, _vbbl_binary_mag_point,
    _vbbl_SG12_5, _adaptive_contouring_linear, _solver)

from MulensModel.pointlens import _AbstractMagnification
from MulensModel.utils import Utils
from MulensModel.version import __version__ as mm_version


class BinaryLens(object):
    """
    *DEPRECATED*

    The BinaryLens class has been replaced with separate classes for each
    magnification method.
    """

    def __init__(self, parameters=None):
        raise NotImplementedError(
            'In version 3, the BinaryLens class has been replaced with a ' +
            'variety of classes with inheritance.')


class _BinaryLensPointSourceMagnification(_AbstractMagnification):
    """
    Equations for calculating point-source--binary-lens magnification.
    This is a placeholder class to establish the basic methods and attributes
    and over-write methods from
    :py:class:`~MulensModel.pointlens.PointSourcePointLensMagnification`
    that do not apply to binary lenses.

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

    """
    def __init__(self, **kwargs):
        super().__init__(trajectory=kwargs['trajectory'])
        self._q = float(self.trajectory.parameters.q)  # This speeds-up code for np.float input.
        self._solver = _solver

        self._source_x = self.trajectory.x
        self._source_y = self.trajectory.y
        self._separations = self.trajectory.parameters.get_s(self.trajectory.times)
        if isinstance(self._separations, (float, int)):
            self._separations = self._separations * np.ones(len(self._source_x))
        self._zip_kwargs = None

    def get_magnification(self):
        """
        Calculate the magnification

        Parameters : None

        Returns :
            magnification: *np.ndarray*
                The magnification for each point in :py:attr:`~trajectory`.
        """
        zip_args = [self._source_x, self._source_y, self._separations]

        out = []
        if self._zip_kwargs is None:
            for (x, y, separation) in zip(*zip_args):
                out.append(self._get_1_magnification(x, y, separation))
        else:
            zip_args += [self._zip_kwargs]
            for (x, y, separation, kwargs_) in zip(*zip_args):
                out.append(self._get_1_magnification(x, y, separation, **kwargs_))

        self._magnification = np.array(out)
        return self._magnification


class _LimbDarkeningForMagnification(object):
    """
    Abstract class for passing information on limb darkening coefficients for magnification calculations
    """
    def _set_LD_coeffs(self, u_limb_darkening, gamma, default_gamma=None):
        """
        Set both u and gamma LD coeffs based on info provided.

        default_gamma should be None or 0.
        """
        if gamma is not None and u_limb_darkening is not None:
            raise ValueError('Only one limb darkening parameter can be set for magnification calculations')
        elif u_limb_darkening is None and gamma is None:
            if default_gamma is None:
                self._gamma = None
                self._u_limb_darkening = None
            else:
                self._gamma = default_gamma
                self._u_limb_darkening = Utils.gamma_to_u(self._gamma)

        elif gamma is None:
            self._u_limb_darkening = float(u_limb_darkening)
            self._gamma = Utils.u_to_gamma(self._u_limb_darkening)
        else:
            self._gamma = float(gamma)
            self._u_limb_darkening = Utils.gamma_to_u(self._gamma)


class BinaryLensPointSourceWM95Magnification(_BinaryLensPointSourceMagnification):
    """
    Equations for calculating point-source--binary-lens magnification following
    the `Witt & Mao 1995, ApJL, 447, L105 <https://ui.adsabs.harvard.edu/abs/1995ApJ...447L.105W/abstract>`_
    prescription assuming a POINT source.

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._mass_1 = 1. / (1. + self._q)
        self._mass_2 = self._q / (1. + self._q)
        self._total_mass = 0.5 * (self._mass_1 + self._mass_2)  # This is total_mass in WM95 paper.
        self._mass_difference = 0.5 * (self._mass_2 - self._mass_1)
        self._position_z2 = 0. + 0.j

        self._position_z1 = None

        self._last_polynomial_input = None
        self._polynomial_roots = None

    def _get_1_magnification(self, x, y, separation):

        """
        Calculate point-source--binary-lens magnification.

        Returns :
            magnification: *float*
                Point source magnification.
        """
        return self._get_1_magnification_point_source(x, y, separation)

    def _get_1_magnification_point_source(self, x, y, separation):
        """
        Calculate point-source--binary-lens magnification.
        """
        (x, y) = self._change_frame(x, y, separation)

        self._zeta = float(x) + float(y) * 1.j
        self._position_z1 = -separation + 0.j

        jacobian_determinant = self._get_jacobian_determinant()
        signed_magnification = 1. / jacobian_determinant
        magnification = fsum(abs(signed_magnification))

        return magnification

    def _change_frame(self, x, y, separation):
        """
        Change frame in which source position is provided:
        from the center of mass to the secondary mass.
        """
        x_new = x - separation / (1. + self._q)
        return (x_new, y)

    def _get_jacobian_determinant(self):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._verify_polynomial_roots())
        add_1 = self._mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self._mass_2 / (self._position_z2 - roots_ok_bar)**2
        derivative = add_1 + add_2
        return 1. - derivative * np.conjugate(derivative)

    def _get_polynomial(self):
        """calculate coefficients of the polynomial in planet frame"""
        total_m = self._total_mass
        m_diff = self._mass_difference
        zeta = self._zeta
        zeta_conj = zeta.conjugate()

        c_sum = Utils.complex_fsum

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

    def _get_polynomial_roots(self):
        """roots of the polynomial"""
        # ***Casting to float speeds-up code for np.float input.***

        polynomial_input = [self._mass_1, self._mass_2, self._position_z1, self._zeta]
        if polynomial_input == self._last_polynomial_input:
            return self._polynomial_roots

        polynomial = self._get_polynomial()

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

    def _verify_polynomial_roots(self, return_distances=False):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots()

        # Two lines below are simplified assuming self._position_z1.imag = 0 and same for z2.
        roots_conj = np.conjugate(roots)
        solutions = (self._zeta +  # This backs-up the lens equation.
                     self._mass_1 / (roots_conj - self._position_z1) + self._mass_2 / (roots_conj - self._position_z2))

        out = []
        distances = []  # The values in distances[] are a diagnostic on how good the numerical accuracy is.
        for (i, root) in enumerate(roots):
            distances_from_root = abs((solutions - root) ** 2)
            min_distance_arg = np.argmin(distances_from_root)

            if i == min_distance_arg:
                out.append(root)
                distances.append(distances_from_root[min_distance_arg])

        # If the lens equation is solved correctly, there should be 3 or 5 solutions (corresponding to 3 or 5 images).
        if len(out) not in [3, 5]:
            self._raise_error_wrong_number_of_solutions(n_solutions=len(out))

        if return_distances:
            return (np.array(out), np.array(distances))
        else:
            return np.array(out)

    def _raise_error_wrong_number_of_solutions(self, n_solutions):
        """Format error message and raise error for wrong number of solutions"""
        separation = -self._position_z1.real
        msg = ("Wrong number of solutions to the lens equation of binary lens.\nGot {:} and expected 3 or 5.\nThe "
               "parameters (m1, m2, s, source_x, source_y, solver) are:\n{:} {:} {:} {:} {:}  {:}\n\n"
               "Consider using 'point_source_point_lens' method for epochs when the source is very far from the lens. "
               "Note that it's different from 'point_source' method.")
        txt = msg.format(n_solutions, repr(self._mass_1), repr(self._mass_2), repr(separation),
                         repr(self._zeta.real), repr(self._zeta.imag), self._solver)

        if self._solver != "Skowron_and_Gould_12":
            txt += ("\n\nYou should switch to using Skowron_and_Gould_12  polynomial root solver. It is much more "
                    "accurate than numpy.polynomial.polynomial.polyroots(). Skowron_and_Gould_12 method is selected "
                    "in automated way if VBBL is imported properly.")

        distance = sqrt(self._zeta.real**2 + self._zeta.imag**2)
        if distance > 200.:
            txt += ("\n\nYou try to calculate magnification for a huge distance between source and lens, "
                    "which is causing the error.")
        elif separation > 100.:
            txt += "\n\nYou try to calculate magnification for a huge separation, which is causing the error."
        elif self._mass_2 > 1.e-6 * self._mass_1 and (distance < 15. or distance < 2. * separation):
            txt += "\n\nThis is a surprising error - please contact code authors and provide the above error message."

        txt += "\nMulensModel version: {:}".format(mm_version)
        raise ValueError(txt)


class BinaryLensPointSourceVBBLMagnification(_BinaryLensPointSourceMagnification):
    """
    Equations for calculating point-source--binary-lens magnification using
    VBBL for point sources.

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`
    """
    def _get_1_magnification(self, x, y, separation):
        """
        Calculate 1 magnification using VBBL.
        """
        return self._get_1_magnification_point_source(float(x), float(y), float(separation))

    def _get_1_magnification_point_source(self, x, y, separation):
        """
        Call VBBL to get 1 magnification for point source.
        This function is also called by child classes.
        """
        return _vbbl_binary_mag_point(separation, self._q, x, y)


class BinaryLensPointSourceMagnification(_BinaryLensPointSourceMagnification):
    """
    Optimal class for calculation of point-source binary-lens magnification:
    it first tries VBBL and then switches to Witt & Mao 1995 if the former fails.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._first = BinaryLensPointSourceVBBLMagnification(**kwargs)
        self._second = BinaryLensPointSourceWM95Magnification(**kwargs)

    def _get_1_magnification(self, x, y, separation):
        """
        Calculate 1 magnification value
        """
        return self._get_1_magnification_point_source(float(x), float(y), float(separation))

    def _get_1_magnification_point_source(self, x, y, separation):
        """
        First try VBBL then WM95 to get magnification
        """
        repeat = False
        try:
            out = self._first._get_1_magnification_point_source(x, y, separation)
        except Exception:
            repeat = True

        if repeat or out < 1.:
            out = self._second._get_1_magnification_point_source(x, y, separation)

        return out


class _FiniteSource(object):
    """
    Abstract class for setting and checking rho value.
    """
    def _set_and_check_rho(self):
        """
        Check if rho is float and positive.
        """
        rho = self.trajectory.parameters.rho
        if rho is None:
            raise TypeError('rho must be positive float, but None was provided')

        if rho < 0:
            raise ValueError('rho must be positive, got: {:}'.format(self._rho))

        self._rho = float(rho)


class BinaryLensQuadrupoleMagnification(BinaryLensPointSourceMagnification, _LimbDarkeningForMagnification,
                                        _FiniteSource):
    """
    Magnification in quadrupole approximation of the binary-lens/finite-source event - based on
    `Gould 2008 ApJ 681, 1593 <https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract>`_.

    The origin of the coordinate system is at the center of mass and
    both masses are on X axis with higher mass at negative X; this
    means that the higher mass is at (X, Y)=(-s*q/(1+q), 0) and
    the lower mass is at (s/(1+q), 0).

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

        gamma: *float*
            Linear limb-darkening coefficient in gamma convention.

    """

    def __init__(self, gamma=None, u_limb_darkening=None, **kwargs):
        super().__init__(**kwargs)
        self._set_LD_coeffs(u_limb_darkening=u_limb_darkening, gamma=gamma, default_gamma=0.)
        self._set_and_check_rho()

        self._point_source_magnification = None
        self._quadupole_magnification = None

    def _get_magnification_w_plus(self, radius, x, y, separation):
        """Evaluates Gould (2008) eq. 7"""
        dx = [1., 0., -1., 0.]
        dy = [0., 1., 0., -1.]
        out = []
        for (dx_, dy_) in zip(dx, dy):
            xx = x + dx_ * radius
            yy = y + dy_ * radius
            out.append(self._get_1_magnification_point_source(x=xx, y=yy, separation=separation))

        return 0.25 * fsum(out) - self._point_source_magnification

    def _get_1_magnification(self, x, y, separation):
        """
        Calculate the magnification

        Parameters : None

        Returns :
            magnification: *float* or *np.ndarray*
                The magnification for each point
                specified by `u` in :py:attr:`~trajectory`.
        """
        # In this function, variables named self._a_* depict magnification.
        separation = float(separation)
        x = float(x)
        y = float(y)
        self._point_source_magnification = self._get_1_magnification_point_source(x=x, y=y, separation=separation)
        self._a_rho_half_plus = self._get_magnification_w_plus(radius=0.5 * self._rho, x=x, y=y, separation=separation)
        self._a_rho_plus = self._get_magnification_w_plus(radius=self._rho, x=x, y=y, separation=separation)

        # This is Gould 2008 eq. 9:
        self._a_2_rho_square = (16. * self._a_rho_half_plus - self._a_rho_plus) / 3.

        # Gould 2008 eq. 6 (part 1/2):
        self._quadrupole_magnification = (
            self._point_source_magnification + self._a_2_rho_square * (1. - 0.2 * self._gamma))

        return self._quadrupole_magnification


class BinaryLensHexadecapoleMagnification(BinaryLensQuadrupoleMagnification):
    """
    Magnification in hexadecapole approximation of the binary-lens/finite-source event - based on
    `Gould 2008 ApJ 681, 1593 <https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract>`_.

    For coordinate system convention see :py:class:`BinaryLensQuadrupoleMagnification`

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

        gamma: *float*
            Linear limb-darkening coefficient in gamma convention.

        all_approximations: *boolean*, optional
            If *True, return hexadecapole, quadrupole, and point source
            approximations. Default is *False*.
    """

    def __init__(self, all_approximations=False, **kwargs):
        super().__init__(**kwargs)
        self._all_approximations = all_approximations

    def _get_1_magnification(self, x, y, separation):
        """
        Calculate the magnification

        Returns :
            magnification: *float* or *sequence* of three *floats*
                Hexadecapole approximation (*float*) by default.
                If :py:attr:`all_approximations` is *True*, returns
                hexadecapole, quadrupole, and point source approximations
                (*sequence* of three *floats*)
        """
        # In this function, variables named a_* depict magnification.
        a_quadrupole = super()._get_1_magnification(x=x, y=y, separation=separation)

        a_rho_times = self._get_magnification_w_times(x=x, y=y, separation=separation)

        # This is Gould (2008) eq. 9:
        a_4_rho_power4 = 0.5 * (self._a_rho_plus + a_rho_times) - self._a_2_rho_square
        # This is Gould (2008) eq. 6 (part 2/2):
        a_add = a_4_rho_power4 * (1. - 11. * self._gamma / 35.)
        a_hexadecapole = a_quadrupole + a_add

        if self._all_approximations:
            return (a_hexadecapole, self._quadrupole_magnification, self._point_source_magnification)
        else:
            return a_hexadecapole

    def _get_magnification_w_times(self, x, y, separation):
        """Evaluates Gould (2008) eq. 8"""
        shift = self._rho / sqrt(2.)
        dx = [1., -1., -1., 1.]
        dy = [1., 1., -1., -1.]
        out = []
        for (dx_, dy_) in zip(dx, dy):
            xx = x + dx_ * shift
            yy = y + dy_ * shift
            out.append(self._get_1_magnification_point_source(x=float(xx), y=float(yy), separation=separation))

        return 0.25 * fsum(out) - self._point_source_magnification


class BinaryLensVBBLMagnification(_BinaryLensPointSourceMagnification, _LimbDarkeningForMagnification, _FiniteSource):
    """
    Binary lens finite source magnification calculated using VBBL
    library that implements advanced contour integration algorithm
    presented by `Bozza 2010 MNRAS, 408, 2188
    <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract>`_.
    See also `VBBL website by Valerio Bozza
    <http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm>`_.

    For coordinate system convention see
    :py:class:`BinaryLensQuadrupoleMagnification`

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

        gamma: *float*, optional
            Linear limb-darkening coefficient in gamma convention.

        u_limb_darkening: *float*, optional
            Linear limb-darkening coefficient in u convention.
            Note that either *gamma* or *u_limb_darkening* can be
            set.  If neither of them is provided then limb
            darkening is ignored.

        accuracy: *float*, optional
            Requested accuracy of the result.

    """

    def __init__(self, gamma=None, u_limb_darkening=None, accuracy=0.001, **kwargs):
        super().__init__(**kwargs)
        self._set_LD_coeffs(u_limb_darkening=u_limb_darkening, gamma=gamma)
        self._set_and_check_rho()

        if accuracy <= 0.:
            raise ValueError(
                "VBBL requires accuracy > 0 e.g. 0.01 or 0.001;" +
                "\n{:} was  provided".format(accuracy))
        self._accuracy = float(accuracy)

        if not _vbbl_wrapped:
            raise ValueError('VBBL was not imported properly')

        if self._u_limb_darkening is None:
            self._vbbl_function = _vbbl_binary_mag_finite
        else:
            self._vbbl_function = _vbbl_binary_mag_dark

    def _get_1_magnification(self, x, y, separation):
        """
        Calculate 1 magnification using VBBL.
        """
        args = [float(separation), self._q, float(x), float(y), self._rho, self._accuracy]
        if self._u_limb_darkening is not None:
            args += [self._u_limb_darkening]

        return self._vbbl_function(*args)


class BinaryLensAdaptiveContouringMagnification(_BinaryLensPointSourceMagnification, _LimbDarkeningForMagnification,
                                                _FiniteSource):
    """
    Binary lens finite source magnification calculated using
    Adaptive Contouring method by `Dominik 2007 MNRAS, 377, 1679
    <https://ui.adsabs.harvard.edu/abs/2007MNRAS.377.1679D/abstract>`_

    See also
    `AdaptiveContouring website by Martin Dominik
    <http://star-www.st-and.ac.uk/~md35/Software.html>`_

    For coordinate system convention see
    :py:func:`point_source_magnification()`

    Arguments :
        trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
            Including trajectory.parameters =
            :py:class:`~MulensModel.modelparameters.ModelParameters`

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

    """

    def __init__(self, gamma=None, u_limb_darkening=None, accuracy=0.1, ld_accuracy=0.001, **kwargs):
        super().__init__(**kwargs)
        self._set_LD_coeffs(u_limb_darkening=u_limb_darkening, gamma=gamma, default_gamma=0.)
        self._set_and_check_rho()

        # Note that this accuracy is not guaranteed.
        if accuracy <= 0.:
            raise ValueError('adaptive_contouring requires accuracy > 0')
        if ld_accuracy <= 0.:
            raise ValueError('adaptive_contouring requires ld_accuracy > 0')

        if not _adaptive_contouring_wrapped:
            raise ValueError('Adaptive Contouring was not imported properly')

        self._accuracy = float(accuracy)
        self._ld_accuracy = float(ld_accuracy)

    def _get_1_magnification(self, x, y, separation):
        """
        Calculate 1 magnification using AC.
        """
        (x, y) = self._change_frame(x, y, separation)
        args = [float(separation), self._q, float(x), float(y),
                self._rho, self._gamma, self._accuracy, self._ld_accuracy]
        return _adaptive_contouring_linear(*args)

    def _change_frame(self, x, y, separation):
        """
        Change frame in which source position is provided:
        mirror in X (i.e., secondary mass has negative X) and mirror in Y.
        """
        return (-x, -y)
