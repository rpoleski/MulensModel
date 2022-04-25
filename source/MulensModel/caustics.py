import numpy as np
from math import cos, sin
import matplotlib.pyplot as plt

from MulensModel.utils import Utils


class Caustics(object):
    """
    Class for the caustic structure corresponding to a given (*q*, *s*),
    i.e. mass ratio and separation. Implemented for 2-body lenses only.

    Attributes :
        q: *float*
            mass ratio between the 2 bodies; always <= 1
        s: *float*
            separation between the 2 bodies (as a fraction of the
            Einstein ring)
    """

    def __init__(self, q, s):
        # Set s, q
        if isinstance(q, (list, np.ndarray)):
            if len(q) > 1:
                raise NotImplementedError(
                    'Too many q. Does not support more than 2 bodies.')
            else:
                q = q[0]
        self.q = q
        self.s = s

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
        return (self._x, self._y)

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
        Solve the caustics polynomial to calculate the critical curve
        and caustic structure.

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
            eiphi = complex(cos(phi), sin(phi))

            # Coefficients of Eq. 6
            coeff_4 = 1.
            coeff_3 = -2. * self.s
            coeff_2 = Utils.complex_fsum([self.s**2, -eiphi])
            coeff_1 = eiphi * (2. * self.s / (1. + self.q))  # The additional
            # parenthesis make it more stable numerically.
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
