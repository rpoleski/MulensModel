import numpy as np
import matplotlib.pyplot as pl

from MulensModel.utils import Utils


class Caustics(object):
    """
    Class for the caustic structure corresponding to a given (q, s),
    i.e. mass ratio and separation.

    Methods:
        get_caustics() - returns x, y vectors for the caustic locations
        plot() - plots the caustic structure using matplotlib.scatter

    Also contains the critical_curve (an internal class).
    """

    def __init__(self, q=None, s=None):
        """
        Create a Caustics object. Both s and q should be set.
        """
        #Set s, q
        self.q = q
        self.s = s

        #Set place holder variables
        self._x = None
        self._y = None
        self._critical_curve = None

    class CriticalCurve(object):
        """
        Internal class of Caustics. Defines the critical curve (in the
        lens plane).
        """

        def __init__(self):
            self.x = []
            self.y = []

    def _calculate(self, n_points=5000):
        """
        Solve the caustics polynomial to calculate the critical curve
        and caustic structure.
        
        Based on Eq. 6 Cassan 2008 modified so origin is center of
        mass and larger mass is on the left. Uses complex coordinates.
        """
        #Initialize variables
        self._x = []
        self._y = []
        self._critical_curve = self.CriticalCurve()

        #Distance between primary mass and center of mass
        xcm_offset = self.q * self.s / (1. + self.q)

        #Solve for the critical curve (and caustic) in complex coordinates.
        #TO DO: Optimize using vectors instead of a loop
        for phi in np.arange(0., 2.*np.pi, 2*np.pi/n_points):
            #Change the angle to a complex number
            x = np.cos(phi)
            y = np.sin(phi)
            eiphi = np.complex(x, y)

            #Coefficients of Eq. 6
            coeff_4 = 1.
            coeff_3 = -2. * self.s
            coeff_2 = Utils.complex_fsum([self.s**2, -eiphi])
            coeff_1 = 2. * self.s * eiphi / (1. + self.q)
            coeff_0 = -self.s**2 * eiphi / (1. + self.q)

            #Find roots
            coeff_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4]
            roots = np.polynomial.polynomial.polyroots(coeff_list)
            #Store results
            for root in roots:
                self._critical_curve.x.append(root.real - xcm_offset)
                self._critical_curve.y.append(root.imag)

                source_plane_position = self._solve_lens_equation(root)
                self._x.append(source_plane_position.real - xcm_offset)
                self._y.append(source_plane_position.imag)

    def get_caustics(self, n_points=5000):
        """
        Returns x and y vectors corresponding to the outlines of the caustics.
        """
        if self._x is None or self._y is None:
            self._calculate(n_points=n_points)
        return self._x, self._y

    def plot(self, n_points=5000, **kwargs):
        """
        Plots the caustics (using matplotlib.scatter()). 
        """
        if self._x is None:
            self._calculate(n_points=n_points)
        pl.scatter(self._x, self._y, **kwargs)

    def _solve_lens_equation(self, complex_value):
        """
        Solve the lens equation for the given point (in complex coordinates).
        """
        complex_conjugate = np.conjugate(complex_value)
        return complex_value - (1. / (1. + self.q)) * (
            (1./complex_conjugate) + (self.q / (complex_conjugate - self.s)) )

    @property
    def critical_curve(self):
        if self._critical_curve is None:
            self._calculate()
        return self._critical_curve

if __name__ == '__main__':
    """
    Testing functions
    """
    import matplotlib.pyplot as pl
    import MulensModel
    q = 0.0053
    s = 0.548

    MODULE_PATH = "/".join(MulensModel.__file__.split("/source")[:-1])
    SAMPLE_FILE_01 = MODULE_PATH + "/data/MB11293_caustics.dat"
    test_caustics = np.genfromtxt(SAMPLE_FILE_01, names=['X', 'Y'], dtype=None)
    offset_x = s * ( -0.5 + (1. / (1. + q)) )
    center_of_mass = s * (q / (1. + q))

    print( 'Center of Mass: ', s * (q / (1. + q)) )

    caustics = Caustics(q=q, s=s)
    x, y = caustics.get_caustics()
    print(x[0:10], y[0:10])

    pl.figure()

    pl.scatter(
        test_caustics['X'], test_caustics['Y'], 
        color='red', label='pyLIMA',alpha=0.3,marker='o')
    caustics.plot(color='blue', label='MulensModel',alpha=0.3,marker='x',
                  n_points=1000.)

    pl.legend(loc='best')

    binary_caustics = Caustics(q=1, s=0.7)
    pl.figure()
    binary_caustics.plot()

    pl.figure()
    pl.scatter(
        binary_caustics.critical_curve.x, binary_caustics.critical_curve.y)
    pl.title('Critical Curves')

    pl.show()

