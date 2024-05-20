import numpy as np
from cmath import sqrt

from MulensModel.causticsbinary import CausticsBinary


class CausticsPointWithShear(CausticsBinary):
    """
    Class for the caustic structure produced by the Chang-Refsdal lens,
    i.e. single mass with shear and convergence.

    Attributes :
        convergence_K: *float*
            convergence of the lens system

        shear_G: *complex*
            shear of the the lens system
    """

    def __init__(self, convergence_K, shear_G):
        # Set K, G
        self.convergence_K = convergence_K
        self.shear_G = shear_G

        # Set place holder variables
        self._x = None
        self._y = None
        self._critical_curve = None

    def _calculate(self, n_points=5000):
        """
        Solve the caustics polynomial to calculate the critical curve
        and caustic structure.
        """
        # Find number of angles so that 4*n_angles is the multiple of 4 that
        # is closest to n_points.
        n_angles = int(n_points / 4. + .5)

        # Initialize variables
        self._x = []
        self._y = []
        self._critical_curve = self.CriticalCurve()

        # Solve for the critical curve (and caustic) in complex coordinates.
        G_conjugate = self.shear_G.conjugate()
        phi = np.linspace(0., 2. * np.pi, n_angles, endpoint=False)
        eiphi = np.exp(1j * phi)
        for eiphi_ in eiphi:
            soln = sqrt(1. / ((1.-self.convergence_K) * eiphi_ + G_conjugate))
            for root in [soln, -soln]:
                self._critical_curve.x.append(root.real)
                self._critical_curve.y.append(root.imag)

                source_plane_position = self._solve_lens_equation(root)
                self._x.append(source_plane_position.real)
                self._y.append(source_plane_position.imag)

    def _solve_lens_equation(self, complex_value):
        """
        Solve the lens equation for the given point (in complex coordinates).
        """
        complex_conjugate = np.conjugate(complex_value)
        return ((1 - self.convergence_K) * complex_value
                - self.shear_G * complex_conjugate - (1. / complex_conjugate))
