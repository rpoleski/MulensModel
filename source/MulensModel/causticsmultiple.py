from warnings import warn
import matplotlib.pyplot as plt
import VBMicrolensing


class CausticMultiple(object):
    """
    Class for the caustic structure corresponding to a given lens geometry,
    i.e. mass ratio and separation. Implemented for 2-body lenses only.

    Attributes :
        geometry: *array* of shape (N x self.n_lenses * 3)
            The geometry of the lens system at given epochs. N = 1
            or N = the number of epochs, if the orbital motion parameters are set.
            for each lens two coridantes and it mass shoule be spefified. For example:
            geometry =    [[0,0,1,            # First lens: x1_1, x1_2, m1
                            1,-0.7,1.e-4,     # Second lens: x2_1, x2_2, m2
                            2,0.7,1.e-4,      # Third lens: x3_1, x3_2, m3
                            0.6,-.6,1.e-6]]   # Fourth lens: x4_1, x4_2, m4
    """

    def __init__(self, geometry):

        self.geometry = geometry
        self._n_lenses = len(geometry[0]) // 3
        if self._n_lenses % 1 != 0:
            raise ValueError("Wrong geometry. Each lens should be specified by 3 numbers: x1, x2, mass.")
        else:
            self._n_lenses = int(self._n_lenses)
        # Set place holder variable
        self._x = None
        self._y = None
        self._critical_curve = None

    def plot(self, n_points=None, plot_lenses=True, **kwargs):
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
        if n_points is not None:
            warn("n_points is not used for multiple lens geometry." +
                 "The number of points is set by VBMicrolensing and cannot be changed.")

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
        if plot_lenses:
            for i in range(self._n_lenses):
                plt.scatter(self.geometry[0][i*3], self.geometry[0][i*3+1], color='k', marker='x')
                plt.annotate(str(i+1), (self.geometry[0][i*3], self.geometry[0][i*3+1]+0.005))

    def get_caustics(self, n_points=None):
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
        if n_points is not None:
            warn("n_points is not used for multiple lens geometry." +
                 "The number of points is set by VBMicrolensing and cannot be changed.")

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
        Calculate the caustic structure for the given lens geometry using VBMicrolensing.
        """
        # Initialize VBMicrolensing() class object
        VBM = VBMicrolensing.VBMicrolensing()
        # Set relative accuracy
        VBM.RelTol = 1e-03
        # Set accuracy
        VBM.Tol = 1e-03
        VBM.SetLensGeometry(self.geometry[0])
        self._critical_curve = self.CriticalCurve()
        caustics = VBM.Multicaustics()
        criticalcurves = VBM.Multicriticalcurves()
        self._x = []
        self._y = []
        for i, one_caustic in enumerate(caustics):
            self._x.extend(one_caustic[0])
            self._y.extend(one_caustic[1])
            self._critical_curve.x.extend(criticalcurves[i][0])
            self._critical_curve.y.extend(criticalcurves[i][1])

    class CriticalCurve(object):
        """
        Internal class of :py:class:`Caustics`. Defines the critical
        curve (in the lens plane). Origin is center of mass of primary and secondary lens components with
        larger mass on the left (*q* < 1).

        Attributes :
            x, y : *list*
                Two lists of length *n_points* giving the x, y
                coordinates of the caustic points.
        """

        def __init__(self):
            self.x = []
            self.y = []
