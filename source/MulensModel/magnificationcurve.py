import math
import warnings

import numpy as np

from MulensModel.binarylens import BinaryLens
from MulensModel.binarylenswithshear import BinaryLensWithShear
from MulensModel.modelparameters import ModelParameters
from MulensModel.pointlens import PointLens, get_pspl_magnification
from MulensModel.trajectory import Trajectory


class MagnificationCurve(object):
    """
    The magnification curve calculated from the model light curve.

    The key function is :py:func:`set_magnification_methods`, which
    specifies the method used to calculate the finite source
    magnification and when to use it..

    Arguments :
        times: iterable of *floats*
            the times at which to generate the magnification curve

        parameters: :py:class:`~MulensModel.modelparameters.ModelParameters`
            specifies the microlensing parameters

        parallax: *dict*, optional
            dictionary specifying what parallax effects should be
            used, e.g., ``{'earth_orbital': True, 'satellite': False,
            'topocentric': False}``

        coords: :py:class:`~MulensModel.coordinates.Coordinates`, optional
            sky coordinates of the event

        satellite_skycoord: *Astropy.coordinates.SkyCoord*, optional
            sky coordinates of the satellite specified by the
            ephemerides file. See
            :py:obj:`MulensModel.mulensdata.MulensData.satellite_skycoord`.

        gamma: *float*, optional
            limb darkening coefficient in gamma convention; defaults to 0

    Attributes :
        trajectory: :py:class:`~MulensModel.Trajectory.trajectory`
            Trajectory used to calculate positions of
            the source that are used to calculate magnification values.
    """

    def __init__(self, times, parameters, parallax=None,
                 coords=None, satellite_skycoord=None, gamma=0.):
        # Set times
        self.times = np.asarray(times)

        # Check for ModelParameters and set.
        if isinstance(parameters, ModelParameters):
            self.parameters = parameters
        else:
            raise ValueError(
                'parameters is a required keyword and must be a ' +
                'ModelParameters object')

        # Calculate the source trajectory (i.e. u(t))
        self.trajectory = Trajectory(
            self.times, parameters=parameters, parallax=parallax,
            coords=coords, satellite_skycoord=satellite_skycoord)

        # Initialize the magnification vector
        self._magnification = None

        # Set methods' variables:
        self._methods_epochs = None
        self._methods_names = []
        self._default_method = None
        self._methods_parameters = None

        self._gamma = gamma

    def set_magnification_methods(self, methods, default_method):
        """
        Sets methods used for magnification calculation.

        For available methods, see:
            :py:func:`get_point_lens_magnification`

            and

            :py:func:`get_binary_lens_magnification`

        Parameters :
            methods: *list*
                List that specifies which methods (*str*) should be
                used when (*float* values for Julian dates). Given
                method will be used for times between the times
                between which it is on the list, e.g.,

                .. code-block:: python

                  methods = [
                      2455746., 'Quadrupole', 2455746.6, 'Hexadecapole',
                      2455746.7, 'VBBL', 2455747., 'Hexadecapole',
                      2455747.15, 'Quadrupole', 2455748.]

            default_method: *str*
                Name of the method to be used for epochs outside the ranges
                specified in *methods*.
        """
        self._default_method = default_method
        if methods is None:
            self._methods_epochs = None
            self._methods_names = []
            return

        if not isinstance(methods, list):
            msg = ('MagnificationCurve.set_magnification_methods() ' +
                   'requires a list as a parameter')
            raise TypeError(msg)
        epochs = methods[0::2]
        names = methods[1::2]

        for epoch in epochs:
            if not isinstance(epoch, (float, int)):
                raise TypeError('Wrong epoch: {:}'.format(epoch))
        for method in names:
            if not isinstance(method, str):
                raise TypeError('Wrong method: {:}'.format(method))
        for (e_beg, e_end) in zip(epochs[::2], epochs[1::2]):
            if e_beg >= e_end:
                msg = ('Incorrect epochs provided: {:} and {:} (first should' +
                       ' be earlier)')
                raise ValueError(msg.format(e_beg, e_end))

        self._methods_epochs = np.array(epochs)
        self._methods_names = names

    def set_magnification_methods_parameters(self, methods_parameters):
        """
        Set additional parameters for magnification calculation methods.

        Parameters :
            methods_parameters: *dict*
                Dictionary that for method names (keys) returns dictionary
                in the form of ``**kwargs`` that are passed to given method,
                e.g., ``{'VBBL': {'accuracy': 0.005}}``.

        """
        self._methods_parameters = methods_parameters

    def get_magnification(self):
        """
        Calculate magnification.

        Returns :
            magnification: *np.ndarray*
                Vector of magnifications.

        """
        if self.parameters.rho is not None:
            self._check_for_finite_source_method()

        if self.parameters.n_lenses == 1:
            magnification = self.get_point_lens_magnification()
        elif (self.parameters.n_lenses == 2):
            magnification = self.get_binary_lens_magnification()
        else:
            raise NotImplementedError(
                "magnification for more than 2 lenses not handled yet")
        self._magnification = magnification
        return self._magnification

    def _check_for_finite_source_method(self):
        """
        check if there is method defined that uses finite source
        calculations and warn if not
        """
        methods = self._methods_names + [self._default_method]
        set_ = set(['point_source', 'point_source_point_lens'])
        if len(set(methods)-set_) == 0:
            warnings.warn('no finite-source method is set', UserWarning)
            return

    def get_point_lens_magnification(self):
        """
        Calculate the Point Lens magnification.

        Allowed magnification methods
        (set by :py:func:`set_magnification_methods()`) :
            ``point_source``:
                standard Paczynski equation for a point source/point lens.

            ``finite_source_uniform_Gould94``:
                Uses the `Gould 1994 ApJ, 421L, 71
                <https://ui.adsabs.harvard.edu/abs/1994ApJ...421L..71G/abstract>`_
                prescription assuming a
                *uniform* (and circular) source. This method interpolates
                pre-computed tables. The relative interpolation
                errors are smaller than 10^-4.

            ``finite_source_uniform_Gould94_direct``:
                Same as ``finite_source_uniform_Gould94``, but calculates
                the underlying functions directly
                (i.e., without interpolation).

            ``finite_source_uniform_WittMao94``:
                Uses the `Witt and Mao 1994 ApJ, 430, 505
                <https://ui.adsabs.harvard.edu/abs/1994ApJ...430..505W/abstract>`_
                method assuming a *uniform* (and circular) source. This method
                interpolates pre-computed tables. The relative interpolation
                errors are smaller than 10^-4.

            ``finite_source_LD_WittMao94``:
                Uses the `Witt and Mao 1994 ApJ, 430, 505`_ method and
                integrates multiple annuli to obtain magnification for
                a circular source *including limb-darkening*. For description
                of integration of multiple annuli see, e.g.,
                `Bozza et al. 2018 MNRAS, 479, 5157
                <https://ui.adsabs.harvard.edu/abs/2018MNRAS.479.5157B/abstract>`_.
                This method interpolates pre-computed tables. The relative
                interpolation errors are smaller than 10^-4.

            ``finite_source_LD_Yoo04``:
                Uses the `Yoo et al. 2004 ApJ, 603, 139
                <https://ui.adsabs.harvard.edu/abs/2004ApJ...603..139Y/abstract>`_
                prescription for
                a circular source *including limb-darkening*. This method
                interpolates pre-computed tables. The relative interpolation
                errors are smaller than 10^-4.

            ``finite_source_LD_Yoo04_direct``:
                Same as ``finite_source_LD_Yoo04``, but calculates
                the underlying functions directly
                (i.e., without interpolation), hence can be slow.

            ``finite_source_uniform_Lee09``:
                Uses the `Lee et al. 2009 ApJ, 695, 200
                <https://ui.adsabs.harvard.edu/abs/2009ApJ...695..200L/abstract>`_
                method for a circular and *uniform* source. This method
                works well for large sources (rho ~ 1).

            ``finite_source_LD_Lee09``:
                Uses the `Lee et al. 2009 ApJ, 695, 200`_ method for
                a circular source *including limb-darkening*. This method
                works well for large sources (rho ~ 1) but can be slow
                compared to other methods.

        Returns :
            magnification: *np.ndarray*
                Vector of magnifications.

        """
        if self.parameters.n_lenses != 1:
            raise ValueError(
                "You're trying to calculate single lens magnification, but "
                "the model provided has " + str(self.parameters.n_lenses) +
                " lenses")

        pspl_magnification = get_pspl_magnification(self.trajectory)
        if self._methods_epochs is None:
            return pspl_magnification
        point_lens = PointLens(self.parameters)
        magnification = pspl_magnification
        u2 = self.trajectory.x**2 + self.trajectory.y**2
        u_all = np.sqrt(u2)
        methods = np.array(self._methods_for_epochs())

        for method in set(methods):
            kwargs = {}
            if self._methods_parameters is not None:
                if method in self._methods_parameters.keys():
                    kwargs = self._methods_parameters[method]
                if kwargs != {}:
                    raise ValueError(
                        'Methods parameters passed, but currently ' +
                        'no point lens method accepts the parameters')
            selection = (methods == method)

            if method.lower() == 'point_source':
                pass  # These cases are already taken care of.
            elif method.lower() == 'finite_source_uniform_Gould94'.lower():
                magnification[selection] = (
                    point_lens.get_point_lens_finite_source_magnification(
                        u=u_all[selection],
                        pspl_magnification=pspl_magnification[selection]))
            elif (method.lower() ==
                  'finite_source_uniform_Gould94_direct'.lower()):
                magnification[selection] = (
                    point_lens.get_point_lens_finite_source_magnification(
                        u=u_all[selection],
                        pspl_magnification=pspl_magnification[selection],
                        direct=True))
            elif method.lower() == 'finite_source_uniform_WittMao94'.lower():
                pl = point_lens
                magnification[selection] = (
                    pl.get_point_lens_large_finite_source_magnification(
                        u=u_all[selection]))
            elif method.lower() == 'finite_source_LD_WittMao94'.lower():
                pl = point_lens
                magnification[selection] = (
                    pl.get_point_lens_large_LD_integrated_magnification(
                        u=u_all[selection], gamma=self._gamma))
            elif method.lower() == 'finite_source_LD_Yoo04'.lower():
                magnification[selection] = (
                    point_lens.get_point_lens_limb_darkening_magnification(
                        u=u_all[selection],
                        pspl_magnification=pspl_magnification[selection],
                        gamma=self._gamma))
            elif method.lower() == 'finite_source_LD_Yoo04_direct'.lower():
                magnification[selection] = (
                    point_lens.get_point_lens_limb_darkening_magnification(
                        u=u_all[selection],
                        pspl_magnification=pspl_magnification[selection],
                        gamma=self._gamma,
                        direct=True))
            elif method.lower() == 'finite_source_uniform_Lee09'.lower():
                magnification[selection] = (
                    point_lens.get_point_lens_uniform_integrated_magnification(
                        u=u_all[selection],
                        rho=self.parameters.rho))
            elif method.lower() == 'finite_source_LD_Lee09'.lower():
                magnification[selection] = (
                    point_lens.get_point_lens_LD_integrated_magnification(
                        u=u_all[selection],
                        rho=self.parameters.rho,
                        gamma=self._gamma))
            else:
                msg = 'Unknown method specified for single lens: {:}'
                raise ValueError(msg.format(method))

        return magnification

    def get_binary_lens_magnification(self):
        """
        Calculate the binary lens magnification.
        If the shear or convergence are set, then they are used.

        Allowed magnification methods
        (set by :py:func:`set_magnification_methods()`) :
            ``point_source``:
                standard point source magnification calculation.

            ``quadrupole``:
                From `Gould 2008 ApJ, 681, 1593
                <https://ui.adsabs.harvard.edu/abs/2008ApJ...681.1593G/abstract>`_.
                See
                :py:func:`~MulensModel.binarylens.BinaryLens.hexadecapole_magnification()`

            ``hexadecapole``:
                From `Gould 2008 ApJ, 681, 1593`_ See
                :py:func:`~MulensModel.binarylens.BinaryLens.hexadecapole_magnification()`

            ``VBBL``:
                Uses VBBinaryLensing (a Stokes theorem/contour
                integration code) by Valerio Bozza
                (`Bozza 2010 MNRAS, 408, 2188
                <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.2188B/abstract>`_).
                See
                :py:func:`~MulensModel.binarylens.BinaryLens.vbbl_magnification()`

            ``Adaptive_Contouring``:
                Uses AdaptiveContouring (a Stokes theorem/contour
                integration code) by Martin Dominik
                (`Dominik 2007 MNRAS, 377, 1679
                <https://ui.adsabs.harvard.edu/abs/2007MNRAS.377.1679D/abstract>`_).
                See
                :py:func:`~MulensModel.binarylens.BinaryLens.adaptive_contouring_magnification()`

                Note that it doesn't work if shear or convergence are set.

            ``point_source_point_lens``:
                Uses point-source _point_-_lens_ approximation; useful when you
                consider binary lens but need magnification very far from
                the lens (e.g. at separation u = 100).

        Returns :
            magnification: *np.ndarray*
                Vector of magnifications.

        """
        if self.parameters.n_lenses != 2:
            raise ValueError(
                "You're trying to calculate binary lens magnification, but "
                "the model provided has " + str(self.parameters.n_lenses) +
                " lenses")

        if not self.parameters.is_external_mass_sheet:
            binary_lens_class = BinaryLens
            kwargs = dict()
        else:
            binary_lens_class = BinaryLensWithShear
            kwargs = {'convergence_K': self.parameters.convergence_K,
                      'shear_G': self.parameters.shear_G}

        out = self._get_binary_lens_magnification(binary_lens_class, kwargs)

        return out

    def _get_binary_lens_magnification(self, binary_lens_class,
                                       optional_kwargs):
        """
        Run binary lens calculation with proper class (binary_lens_class) and
        some kwargs (optional_kwargs of type *dict*).
        """
        q = self.parameters.q
        binary_kwargs = optional_kwargs
        binary_kwargs['mass_1'] = 1. / (1. + q)
        binary_kwargs['mass_2'] = q / (1. + q)

        is_static = self.parameters.is_static()
        if is_static:
            binary_lens = binary_lens_class(separation=self.parameters.s,
                                            **binary_kwargs)
        methods = self._methods_for_epochs()

        magnification = []
        for index in range(len(self.times)):
            x = self.trajectory.x[index]
            y = self.trajectory.y[index]
            method = methods[index].lower()
            if not is_static:
                binary_lens = binary_lens_class(
                    separation=self.parameters.get_s(self.times[index]),
                    **binary_kwargs)

            kwargs = {}
            if self._methods_parameters is not None:
                if method in self._methods_parameters.keys():
                    kwargs = self._methods_parameters[method]
                    if method not in ['vbbl', 'adaptive_contouring']:
                        msg = ('Methods parameters passed for method {:}' +
                               ' which does not accept any parameters')
                        raise ValueError(msg.format(method))

            if method == 'point_source':
                try:
                    m = binary_lens.point_source_magnification(x, y)
                except Exception as e:
                    text = "Model parameters for above exception:\n"
                    text += str(self.parameters)
                    raise ValueError(text) from e
                    # The code above is based on
                    # https://stackoverflow.com/questions/6062576/
                    # adding-information-to-an-exception/6062799
            elif method == 'quadrupole':
                m = binary_lens.hexadecapole_magnification(
                    x, y, rho=self.parameters.rho, quadrupole=True,
                    gamma=self._gamma)
            elif method == 'hexadecapole':
                m = binary_lens.hexadecapole_magnification(
                    x, y, rho=self.parameters.rho, gamma=self._gamma)
            elif method == 'vbbl':
                if isinstance(binary_lens, BinaryLensWithShear):
                    raise ValueError("Finite source VBBL is not available "
                                     "for BinaryLensWithShear")
                m = binary_lens.vbbl_magnification(
                    x, y, rho=self.parameters.rho, gamma=self._gamma, **kwargs)
            elif method == 'adaptive_contouring':
                if isinstance(binary_lens, BinaryLensWithShear):
                    raise ValueError("Adaptive contouring is not available "
                                     "for BinaryLensWithShear")
                m = binary_lens.adaptive_contouring_magnification(
                    x, y, rho=self.parameters.rho, gamma=self._gamma, **kwargs)
            elif method == 'point_source_point_lens':
                u = math.sqrt(x**2 + y**2)
                m = get_pspl_magnification(u)
            else:
                msg = 'Unknown method specified for binary lens: {:}'
                raise ValueError(msg.format(method))

            magnification.append(m)

        return np.array(magnification)

    def _methods_for_epochs(self):
        """
        for given epochs, decide which methods should be used to
        calculate magnification, but don't run the calculations
        """
        out = [self._default_method] * len(self.times)
        if self._methods_epochs is None:
            return out

        brackets = np.searchsorted(self._methods_epochs, self.times)
        n_max = len(self._methods_epochs)

        out = [self._methods_names[value - 1]
               if (value > 0 and value < n_max) else self._default_method
               for value in brackets]
        return out
