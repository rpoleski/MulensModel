import warnings
import numpy as np
from math import fsum

from MulensModel.utils import Utils


class TripleLens(object):
    """
    The triple lens equation - its solutions, images, parities,
    magnifications, etc.

    Arguments :
        mass_1: *float*
            Mass of the primary as a fraction of the total mass.

        mass_2: *float*
            Mass of the secondary as a fraction of the total mass.

        mass_3: *float*
            Mass of the tertiary as a fraction of the total mass.
            XXX NOTE - add conventions

        separation_21: *float*
            Separation between bodies 2 and 1 as a fraction of
            the Einstein ring of the whole system.

        separation_31: *float*
            Separation between bodies 2 and 1 as a fraction of
            the Einstein ring of the whole system.

# DECIDE IF IT'S psi OR phi AND APPLY TO WHOLE CODE
        psi: *float*
            XXX NOTE - description

    Note: masses 1, 2, and 3 may be defined as a fraction of some other
    mass than the total mass. This is possible but not recommended -
    make sure you know what you're doing before you start using this
    possibility.
    """
    def __init__(self, mass_1=None, mass_2=None, mass_3=None,
                 separation_21=None, separation_31=None, psi=None):
        self._mass_1 = mass_1
        self._mass_2 = mass_2
        self._mass_3 = mass_3
        # XXX warning if sum doesn't ~1
        self._separation_21 = separation_21
        self._separation_31 = separation_31
        self._psi = psi

        positions = Utils._parameters_to_center_of_mass_coords_3L(
            q_21=self._mass_2/self._mass_1, q_31=self._mass_3/self._mass_1,
            s_21=self._separation_21, s_31=self._separation_31, psi=self._psi)
        self._position_z1 = positions[0, 0] + positions[0, 1] * 1.j
        self._position_z2 = positions[1, 0] + positions[1, 1] * 1.j
        self._position_z3 = positions[2, 0] + positions[2, 1] * 1.j

        self._solver = 'numpy'  # XXX

    def get_point_source_magnification(self, source_x, source_y):
        """
        XXX

        Parameters :
            source_x: *float*
                X-axis coordinate of the source.

            source_y: *float*
                Y-axis coordinate of the source.

        Returns :
            magnification: *float*
                Point-source triple-lens magnification.
        """
        # use Rhie 2002 https://arxiv.org/abs/astro-ph/0202294
        signed_magnification = 1. / self._jacobian_determinant_ok(
            source_x=source_x, source_y=source_y)
        return fsum(abs(signed_magnification))

    def _get_polynomial_roots(self, source_x, source_y):
        """roots of the polynomial"""
        # XXX - lazy loading like in binarylens.py?

        polynomial = self._get_R02_polynomial(
            source_x=source_x, source_y=source_y)

        np_polyroots = np.polynomial.polynomial.polyroots
        if self._solver == 'numpy':
            polynomial_roots = np_polyroots(polynomial)
        elif self._solver == 'Skowron_and_Gould_12':
            raise NotImplementedError('We have not yet implemented loading ' +
                                      'SG12 from binarylens.py')  # XXX
        else:
            raise ValueError('Unknown solver: {:}'.format(self._solver))

        return polynomial_roots

    def _polynomial_roots_ok(self, source_x, source_y):
        """verified roots of polynomial i.e. roots of lens equation"""
        roots = self._get_polynomial_roots(
            source_x=source_x, source_y=source_y)

        add_1 = self._mass_1 / np.conjugate(roots - self._position_z1)
        add_2 = self._mass_2 / np.conjugate(roots - self._position_z2)
        add_3 = self._mass_3 / np.conjugate(roots - self._position_z3)
        solutions = source_x + source_y * 1.j + add_1 + add_2 + add_3

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

        if len(out) not in [4, 6, 8, 10]:
            warnings.warn(str(len(out)) +  # XXX - warning here
                          " solustions for triple lens equation", UserWarning)

        return np.array(out)

    def _jacobian_determinant_ok(self, source_x, source_y):
        """determinants of lens equation Jacobian for verified roots"""
        roots_ok_bar = np.conjugate(self._polynomial_roots_ok(
                                    source_x=source_x, source_y=source_y))
        # Variable X_bar is conjugate of variable X.
        add_1 = self._mass_1 / (self._position_z1 - roots_ok_bar)**2
        add_2 = self._mass_2 / (self._position_z2 - roots_ok_bar)**2
        add_3 = self._mass_3 / (self._position_z3 - roots_ok_bar)**2
        derivative = add_1 + add_2 + add_3

        return 1. - derivative * np.conjugate(derivative)

    def _get_R02_polynomial(self, source_x, source_y):
        """
        Calculate polynomial coefficients using Rhie (2002).
        """
        omega = source_x + source_y * 1.j

        x_1 = self._position_z1
        x_2 = self._position_z2
        x_3 = self._position_z3

        epsilon_1 = self._mass_1
        epsilon_2 = self._mass_2
        epsilon_3 = self._mass_3

        a = -x_1 - x_2 - x_3
        b = x_1 * x_2 + x_1 * x_3 + x_2 * x_3
        c = -x_1 * x_2 * x_3
        d = (epsilon_1 * x_2 * x_3 + epsilon_2 * x_1 * x_3 +
             epsilon_3 * x_1 * x_2)

        H_0 = np.zeros(11, dtype=np.complex_)
        H_1 = np.zeros(11, dtype=np.complex_)
        H_2 = np.zeros(11, dtype=np.complex_)
        H_3 = np.zeros(11, dtype=np.complex_)
        # We define these up to H_x[10] so that equation for cff[10] is
        # the same as for other indexes.

        H_3[9] = 1
        H_3[8] = 3 * a
        H_3[7] = 3 * b + 3 * a**2
        H_3[6] = 3 * c + 6 * a * b + a**3
        H_3[5] = 6 * a * c + 3 * b**2 + 3 * a**2 * b
        H_3[4] = 6 * b * c + 3 * a**2 * c + 3 * a * b**2
        H_3[3] = 3 * c**2 + 6 * a * b * c + b**3
        H_3[2] = 3 * a * c**2 + 3 * b**2 * c
        H_3[1] = 3 * b * c**2
        H_3[0] = c**3
        H_2[8] = 1
        H_2[7] = 3 * a
        H_2[6] = d + 2 * b + 3 * a**2
        H_2[5] = 2 * a * d + 4 * a * b + a**3 + 2 * c
        H_2[4] = 2 * d * b + d * a**2 + 4 * a * c + 2 * a**2 * b + b**2
        H_2[3] = (
            2 * d * c + 2 * d * a * b + 2 * a**2 * c + a * b**2 + 2 * b * c)
        H_2[2] = 2 * c * a * d + d * b**2 + 2 * a * b * c + c**2
        H_2[1] = 2 * b * c * d + a * c**2
        H_2[0] = c**2 * d
        H_1[7] = 1
        H_1[6] = 3 * a
        H_1[5] = 2 * d + 3 * a**2 + b
        H_1[4] = 4 * a * d + a**3 + 2 * a * b + c
        H_1[3] = d**2 + 2 * a**2 * d + 2 * b * d + b * a**2 + 2 * a * c
        H_1[2] = a * d**2 + 2 * a * b * d + 2 * c * d + c * a**2
        H_1[1] = b * d**2 + 2 * a * c * d
        H_1[0] = c * d**2
        H_0[6] = 1
        H_0[5] = 3 * a
        H_0[4] = 3 * d + 3 * a**2
        H_0[3] = 6 * a * d + a**3
        H_0[2] = 3 * d**2 + 3 * a**2 * d**2
        H_0[1] = 3 * a * d**2
        H_0[0] = d**3

        omega_bar = np.conjugate(omega)
        omega_1_bar = omega_bar - np.conjugate(x_1)
        omega_2_bar = omega_bar - np.conjugate(x_2)
        omega_3_bar = omega_bar - np.conjugate(x_3)
        a_omega = omega_1_bar + omega_2_bar + omega_3_bar
        b_omega = (omega_1_bar * omega_2_bar +
                   omega_2_bar * omega_3_bar +
                   omega_3_bar * omega_1_bar)
        c_omega = omega_1_bar * omega_2_bar * omega_3_bar

        cff = np.zeros(11, dtype=np.complex_)
        # XXX polynomial convention: cff[i] * z**i
        for k in range(1, 10+1):
            cff[k] = (
                H_0[k-1] + H_1[k-1] * a_omega + H_2[k-1] * b_omega +
                H_3[k-1] * c_omega -
                H_0[k] * omega - H_1[k] * (omega * a_omega - 1) -
                H_2[k] * (omega * b_omega + a_omega - omega_bar) -
                H_3[k] * (omega * c_omega + b_omega))
        cff[0] = (  # XXX - it's not totally clear if that what Eq. 7 means:
            -H_0[0] * omega - H_1[0] * (omega * a_omega - 1) +
            -H_2[0] * (omega * b_omega + a_omega - omega_bar) +
            -H_3[0] * (omega * c_omega + b_omega))

        return cff

# XXX :
#    def get_hexadecapole_magnification(self, source_x, source_y, rho, gamma,
#                               quadrupole=False, all_approximations=False):
