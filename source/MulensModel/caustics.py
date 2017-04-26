import numpy as np

from MulensModel import Utils

### NEEDS COMMENTS!!!!
class Caustics(object):

    def __init__(self, q=None, s=None):
        self.q = q
        self.s = -s
        self._x = None
        self._y = None

    def calculate(self, npoints=5000):
        """
        coefficients from Eq. 6 Cassan 2008. Larger mass at the origin.
        # Change so default origin is center of mass
        """
        self._x = []
        self._y = []
        for phi in np.arange(0., 2.*np.pi, 2*np.pi/npoints):
            x = np.cos(phi)
            y = np.sin(phi)
            eiphi = np.complex(x, y)

            coeff_4 = 1.
            coeff_3 = 2. * self.s
            coeff_2 = Utils.complex_fsum([self.s**2, -eiphi])
            coeff_1 = -2. * self.s * eiphi / (1. + self.q)
            coeff_0 = -self.s**2 * eiphi / (1. + self.q)

            coeff_list = [coeff_0, coeff_1, coeff_2, coeff_3, coeff_4]
            roots = np.polynomial.polynomial.polyroots(coeff_list)

            for root in roots:
                source_plane_position = self.solve_lens_equation(root)
                self._x.append(source_plane_position.real)
                self._y.append(source_plane_position.imag)

    def plot(self, npoints=5000, **kwargs):
        if self._x is None:
            self.calculate(npoints=npoints)
        pl.scatter(self._x, self._y, **kwargs)

    def solve_lens_equation(self, complex_value):
        complex_conjugate = np.conjugate(complex_value)
        return complex_value - (1. / (1. + self.q)) * (
            (1./complex_conjugate) + (self.q / (complex_conjugate + self.s)) )

if __name__ == '__main__':
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

    pl.figure()

    pl.scatter(
        test_caustics['X']+offset_x+center_of_mass, test_caustics['Y'], 
        color='red', label='Fortran')
    caustics.plot(color='blue', label='MulensModel')

    pl.legend(loc='best')
    pl.show()
