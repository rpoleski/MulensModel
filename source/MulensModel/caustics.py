import numpy as np
from math import cos, sin
import matplotlib.pyplot as plt

from MulensModel.utils import Utils


class Caustics(object):
    """
    Class for the caustic structure corresponding to parameters of either
    binary lens (*q*, *s*) or
    triple lens (*q_21*, *q_31*, *s_21*, *s_31*, *psi*).

    Arguments :
        q: *float*
            mass ratio between the 2 bodies; always <= 1

        s: *float*
            separation between the 2 bodies (as a fraction of the
            Einstein ring)

        q_21: *float*
            mass 2 divided by mass 1

        q_31: *float*
            mass 3 divided by mass 1

        s_21: *float*
            separation between body 2 and body 1

        s_31: *float*
            separation between body 3 and body 1

        psi: *float*
            angle between *s_21* and *s_31*
    """

    def __init__(self,
                 q=None, s=None,
                 q_21=None, q_31=None, s_21=None, s_31=None, psi=None):
        key = {
            'q': q is not None, 's': s is not None, 'q_21': q_21 is not None,
            'q_31': q_31 is not None, 's_21': s_21 is not None,
            's_31': s_31 is not None, 'psi': psi is not None}
        keys_3L = ['q_21', 'q_31', 's_21', 's_31', 'psi']
        if key['q'] and key['s']:
            if np.any([key[k] for k in keys_3L]):
                raise ValueError(
                    'If you provide binary lens parameters (s, q), then do ' +
                    'not provide triple lens parameters: ' +
                    str([k for k in key if key[k] and k in keys_3L]))
            self._n_lenses = 2
            self.q = q  # XXX NEEDS DOCSTRING?
            self.s = s
        elif np.all([key[k] for k in keys_3L]):
            if key['q'] or key['s']:
                raise ValueError(
                    'If you provide triple lens parameters ' +
                    '(q_21, q_31, s_21, s_31, psi) then do not provide ' +
                    'binary lens parameters (s and q).')
            self._n_lenses = 3
            self._q_21 = q_21
            self._q_31 = q_31
            self._s_21 = s_21
            self._s_31 = s_31
            self._psi = psi
        else:
            raise ValueError(
                'Not enough information to define binary or triple lens.' +
                'Provided parameters: ' +
                str([k for k in key if key[k] and k]))

        # Set place holder variables
        self._x = None
        self._y = None
        self._critical_curve = None

    def plot(self, n_points=5000, **kwargs):
        """
        Plots the caustics using :py:func:`matplotlib.pyplot.scatter()`.

        Parameters :
            n_points: *int*, optional
                The number of points to calculate along the caustic.
                Defaults to 5000.

            ``**kwargs``:
                keywords accepted by :py:func:`matplotlib.pyplot.scatter()`

        Note that default scaling of axis may not be equal on both axis.
        To mitigate this, use:
        ``plt.gca().set_aspect('equal')`` or ``plt.axis('equal')``
        (the other possible options are ``'scaled'`` or ``'square'``).

        """
        if "linewidths" not in kwargs and "lw" not in kwargs:
            kwargs["lw"] = 0.
        if self._x is None or len(self._x) != n_points:
            self._calculate(n_points=n_points)

        try:
            plt.scatter(self._x, self._y, **kwargs)
        except Exception:
            print("kwargs passed to plt.scatter():")
            print(kwargs)
            raise

    def get_caustics(self, n_points=5000):
        """
        Returns x and y vectors corresponding to the outlines of the
        caustics.  Origin is center of mass and larger mass is on the
        left (for *q* < 1).

        Parameters:
            n_points : *int*, optional
                The number of points to calculate along the caustic.

        Returns:
            x, y : *list*
                Two lists of length *n_points* giving the *x*, *y*
                coordinates of the caustic points.
        """
        if self._x is None or self._y is None:
            self._calculate(n_points=n_points)
        return self._x, self._y

    @property
    def critical_curve(self):
        """
        Critical curve stored as :py:class:`CriticalCurve` object, read-only
        """
        if self._critical_curve is None:
            self._calculate()
        return self._critical_curve

    def _calculate(self, n_points=5000):
        """
        calculate critical curve and caustic curve
        """
        if self._n_lenses == 2:
            self._calculate_binary_lens(n_points)
        elif self._n_lenses == 3:
            self._calculate_triple_lens(n_points)
        else:
            raise ValueError(
                'Wrong number of lenses: {:}'.format(self._n_lenses))

    def _calculate_binary_lens(self, n_points):
        """
        Solve the caustics polynomial to calculate the critical curve
        and caustic structure for triple lens.

        Based on Eq. 6 Cassan 2008 modified so origin is center of
        mass and larger mass is on the left. Uses complex coordinates.
        """
        # Find number of angles so that 4*n_angles is the multiple of 4 that
        # is closest to n_points.
        n_angles = int(n_points/4.+.5)

        # Initialize variables
        self._x = []
        self._y = []
        self._critical_curve = self.CriticalCurve()

        # Distance between primary mass and center of mass
        xcm_offset = self.q * self.s / (1. + self.q)

        # Solve for the critical curve (and caustic) in complex coordinates.
        for phi in np.linspace(0., 2.*np.pi, n_angles, endpoint=False):
            # Change the angle to a complex number
            eiphi = np.complex(cos(phi), sin(phi))

            # Coefficients of Eq. 6
            coeff_4 = 1.
            coeff_3 = -2. * self.s
            coeff_2 = Utils.complex_fsum([self.s**2, -eiphi])
            coeff_1 = 2. * self.s * eiphi / (1. + self.q)
            coeff_0 = -self.s**2 * eiphi / (1. + self.q)

            # Find roots
            coeff_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4]
            roots = np.polynomial.polynomial.polyroots(coeff_list)
            # Store results
            for root in roots:
                self._critical_curve.x.append(root.real - xcm_offset)
                self._critical_curve.y.append(root.imag)

                source_plane_position = self._solve_lens_equation(root)
                self._x.append(source_plane_position.real - xcm_offset)
                self._y.append(source_plane_position.imag)

    def _solve_lens_equation(self, complex_value):
        """
        Solve the lens equation for the given point (in complex coordinates).
        """
        complex_conjugate = np.conjugate(complex_value)
        return complex_value - (1. / (1. + self.q)) * (
            (1./complex_conjugate) + (self.q / (complex_conjugate - self.s)))

    def _calculate_triple_lens(self, n_points):
        """
        Solve the caustics polynomial to calculate the critical curve
        and caustic structure for triple lens.
        """
        xy = Utils._parameters_to_center_of_mass_coords_3L(
            self._q_21, self._q_31, self._s_21, self._s_31, self._psi)
        z_1 = np.complex(xy[0, 0], x[0, 1])
        z_2 = np.complex(xy[1, 0], x[1, 1])
        z_3 = np.complex(xy[2, 0], x[2, 1])

        epsilon_1 = 1. / (1. + self._q_21 + self._q_31)
        epsilon_2 = epsilon_1 * self._q_21
        epsilon_3 = epsilon_1 * self._q_31

        self._calculate_triple_lens_complex(
            z_1, z_2, z_3, epsilon_1, epsilon_2, epsilon_3, n_points)

    def _calculate_triple_lens_complex(
            self, z_1, z_2, z_3, epsilon_1, epsilon_2, epsilon_3, n_points):
        """
        Calculate triple lens caustic and critical curve using complex
        coordinates as input.
        """
        n_angles = int(n_points/6.+.5)

        # Initialize variables
        self._x = []
        self._y = []
        self._critical_curve = self.CriticalCurve()

        z_1_bar = np.conjugate(z_1)
        z_2_bar = np.conjugate(z_2)
        z_3_bar = np.conjugate(z_3)

        for phi in np.linspace(0., 2.*np.pi, n_angles, endpoint=False):
            polynomial = self._polynomial_3L(
                z_1, z_2, z_3, epsilon_1, epsilon_2, epsilon_3, self._phi)
            roots = np.polynomial.polynomial.polyroots(polynomial)
            for root in roots:
                self._critical_curve.x.append(root.real)
                self._critical_curve.y.append(root.imag)
                root_bar = np.conjugate(root)
                z = (root + epsilon_1 / (z_1_bar - root_bar) +
                     epsilon_2 / (z_2_bar - root_bar) +
                     epsilon_3 / (z_3_bar - root_bar))
                self._x.append(z.real)
                self._y.append(z.imag)

    class CriticalCurve(object):
        """
        Internal class of :py:class:`Caustics`. Defines the critical
        curve (in the lens plane). Origin is center of mass with
        larger mass on the left (*q* < 1).

        Attributes :
            x, y : *list*
                Two lists of length *n_points* giving the x, y
                coordinates of the caustic points.

        """

        def __init__(self):
            self.x = []
            self.y = []
