import warnings
from os.path import join
import numpy as np

import MulensModel as mm


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
        self.times = np.atleast_1d(times)

        # Check for ModelParameters and set.
        if isinstance(parameters, mm.ModelParameters):
            self.parameters = parameters
        else:
            raise TypeError(
                'parameters is a required keyword and must be a ' +
                'ModelParameters object')

        # Trajectory parameters:
        self.parallax = parallax
        self.coords = coords
        self.satellite_skycoord = satellite_skycoord

        # Initialize the magnification vector
        self._magnification = None
        self._magnification_objects = None

        # Set methods' variables:
        self._methods_epochs = None
        self._methods_names = []
        self._default_method = 'point_source'
        self._methods_parameters = None
        self._methods_for_epochs = None
        self._methods_indices = None

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

        epochs, names = self._raise_errors_set_magnification_methods(methods)
        self._methods_epochs = np.array(epochs)
        self._methods_names = names

        self._methods_for_epochs = None
        self._methods_indices = None
        _ = self.methods_for_epochs
        _ = self.methods_indices

    def _raise_errors_set_magnification_methods(self, methods):
        """
        Check if the list of methods is valid and raise errors if not.
        """
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

        return epochs, names

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
        if self.parameters.n_lenses == 1:
            magnification = self.get_point_lens_magnification()
        elif self.parameters.n_lenses == 2:
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
        set_ = set(['point_source', 'point_source_point_lens', None])
        if len(set(methods)-set_) == 0:
            path = join(
                mm.MODULE_PATH, "documents", "magnification_methods.pdf")
            msg = ("A finite source parameter (rho or t_star) is set and no "
                   "finite-source method is set.\n"
                   "For possible magnification methods see\n" +
                   path + "\nor\n"
                   "https://github.com/rpoleski/MulensModel/blob/master/"
                   "documents/magnification_methods.pdf")
            warnings.warn(msg, UserWarning)
            return

    def _set_magnification_objects(self):
        """ High-level function that separations PL/BL and shear/no shear """
        if self.parameters.n_lenses == 1:
            if not self.parameters.is_external_mass_sheet:
                self._set_point_lens_magnification_objects()
            else:
                self._set_point_lens_w_shear_magnification_objects()
        elif self.parameters.n_lenses == 2:
            if not self.parameters.is_external_mass_sheet:
                self._set_binary_lens_magnification_objects()
            else:
                self._set_binary_lens_w_shear_magnification_objects()

    def _setup_trajectory(self, selection):
        """ Create a trajectory object for a given subset of the data
        specified by *selection*. """
        if self.satellite_skycoord is not None:
            satellite_skycoord = self.satellite_skycoord[selection]
        else:
            satellite_skycoord = None

        trajectory = mm.Trajectory(
            self.times[selection], parameters=self.parameters,
            parallax=self.parallax, coords=self.coords,
            satellite_skycoord=satellite_skycoord)
        return trajectory

    def _setup_kwargs(self, method):
        """ Setup the kwargs for a given magnification object."""
        kwargs = {}
        if self._methods_parameters is not None:
            if method.lower() in self._methods_parameters.keys():
                kwargs = self._methods_parameters[method.lower()]

        return kwargs

    def _set_point_lens_magnification_objects(self):
        """ For simple point lens models, create a *dict* of magnification
        objects corresponding to the user-specified magnification methods."""
        self._magnification_objects = {}
        for method, selection in self.methods_indices.items():
            trajectory = self._setup_trajectory(selection)
            kwargs = self._setup_kwargs(method)

            if kwargs != {}:
                raise ValueError(
                    'Methods parameters passed, but currently ' +
                    'no point lens method accepts the parameters')

            if method.lower() == 'point_source':
                self._magnification_objects[method] = \
                    mm.pointlens.PointSourcePointLensMagnification(
                        trajectory=trajectory)
            elif method.lower() == 'finite_source_uniform_Gould94'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceUniformGould94Magnification(
                        trajectory=trajectory)
            elif (method.lower() ==
                  'finite_source_uniform_Gould94_direct'.lower()):
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceUniformGould94Magnification(
                        trajectory=trajectory, direct=True)
            elif method.lower() == 'finite_source_LD_Yoo04'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceLDYoo04Magnification(
                        trajectory=trajectory, gamma=self._gamma)
            elif method.lower() == 'finite_source_LD_Yoo04_direct'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceLDYoo04Magnification(
                        trajectory=trajectory, gamma=self._gamma,
                        direct=True)
            elif method.lower() == 'finite_source_uniform_WittMao94'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceUniformWittMao94Magnification(
                    trajectory=trajectory)
            elif method.lower() == 'finite_source_LD_WittMao94'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceLDWittMao94Magnification(
                    trajectory=trajectory, gamma=self._gamma)
            elif method.lower() == 'finite_source_uniform_Lee09'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceUniformLee09Magnification(
                        trajectory=trajectory)
            elif method.lower() == 'finite_source_LD_Lee09'.lower():
                self._magnification_objects[method] = \
                    mm.pointlens.FiniteSourceLDLee09Magnification(
                        trajectory=trajectory, gamma=self._gamma)
            else:
                msg = 'Unknown method specified for single lens: {:}'
                raise ValueError(msg.format(method))

    def _set_point_lens_w_shear_magnification_objects(self):
        """ For point lens + shear models, create a *dict* of magnification
        objects corresponding to the user-specified magnification methods."""
        self._magnification_objects = {}
        for method, selection in self.methods_indices.items():
            trajectory = self._setup_trajectory(selection)
            kwargs = self._setup_kwargs(method)

            if kwargs != {}:
                raise ValueError(
                    'Methods parameters passed, but currently ' +
                    'no point lens method accepts the parameters')

            if method.lower() == 'point_source':
                self._magnification_objects[method] = \
                    mm.PointSourcePointLensWithShearMagnification(
                        trajectory=trajectory)
            else:
                msg = 'Unknown method specified for single lens with shear: {:}'
                raise ValueError(msg.format(method))

    def get_point_lens_magnification(self):
        """
        Calculate the Point Lens magnification.

        Allowed magnification methods (set by :py:func:`set_magnification_methods()`) :
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

        return self._get_magnification_universal()

    def _set_binary_lens_magnification_objects(self):
        """ For simple binary lens models, create a *dict* of magnification
        objects corresponding to the user-specified magnification methods."""
        self._magnification_objects = {}
        for method, selection in self.methods_indices.items():
            trajectory = self._setup_trajectory(selection)
            kwargs = self._setup_kwargs(method)

            if ((kwargs != {}) and
                    (method.lower() not in ['vbbl', 'adaptive_contouring'])):
                msg = ('Methods parameters passed for method {:}' +
                       ' which does not accept any parameters')
                raise ValueError(msg.format(method))

            if method.lower() == 'point_source':
                self._magnification_objects[method] = \
                    mm.binarylens.BinaryLensPointSourceMagnification(trajectory=trajectory)
            elif method.lower() == 'quadrupole':
                self._magnification_objects[method] = \
                    mm.binarylens.BinaryLensQuadrupoleMagnification(
                    trajectory=trajectory, gamma=self._gamma)
            elif method.lower() == 'hexadecapole':
                self._magnification_objects[method] = \
                    mm.binarylens.\
                    BinaryLensHexadecapoleMagnification(
                    trajectory=trajectory, gamma=self._gamma)
            elif method.lower() == 'vbbl':
                self._magnification_objects[method] = \
                    mm.binarylens. \
                    BinaryLensVBBLMagnification(
                        trajectory=trajectory, gamma=self._gamma, **kwargs)
            elif method.lower() == 'adaptive_contouring':
                self._magnification_objects[method] = \
                    mm.binarylens. \
                    BinaryLensAdaptiveContouringMagnification(
                        trajectory=trajectory, gamma=self._gamma, **kwargs)
            elif method.lower() == 'point_source_point_lens':
                self._magnification_objects[method] = \
                    mm.pointlens.PointSourcePointLensMagnification(
                        trajectory=trajectory)
            else:
                msg = 'Unknown method specified for binary lens: {:}'
                raise ValueError(msg.format(method))

    def _set_binary_lens_w_shear_magnification_objects(self):
        """
        For binary lens + shear models, create a *dict* of magnification
        objects corresponding to the user-specified magnification methods.
        """
        self._magnification_objects = {}
        for method, selection in self.methods_indices.items():
            trajectory = self._setup_trajectory(selection)
            K = self.parameters.parameters.get('convergence_K', 0)
            G = self.parameters.parameters.get('shear_G', complex(0, 0))
            kwargs = {'convergence_K': K, 'shear_G': G}

            if method.lower() == 'point_source':
                self._magnification_objects[method] = \
                    mm.binarylenswithshear. \
                    BinaryLensPointSourceWithShearVBBLMagnification(
                            trajectory=trajectory, **kwargs)
            elif method.lower() == 'point_source_wm95':
                self._magnification_objects[method] = \
                    mm.binarylenswithshear. \
                    BinaryLensPointSourceWithShearWM95Magnification(
                            trajectory=trajectory, **kwargs)
            else:
                msg = 'Unknown method specified for binary lens: {:}'
                raise ValueError(msg.format(method))

    def get_binary_lens_magnification(self):
        """
        Calculate the binary lens magnification.
        If the shear or convergence are set, then they are used.

        Allowed magnification methods (set by :py:func:`set_magnification_methods()`) :
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

        return self._get_magnification_universal()

    def _get_magnification_universal(self):
        """
        Calculate the magnification both for point and binary lens models,
        avoiding code duplication.
        """
        if self.parameters.is_finite_source():
            self._check_for_finite_source_method()

        if self._magnification_objects is None:
            self._set_magnification_objects()

        magnification = np.zeros(len(self.times))
        for method, selection in self.methods_indices.items():
            magnification[selection] = \
                self._magnification_objects[method].get_magnification()

        return magnification

    def get_d_A_d_params(self, parameters):
        """
        Calculate d A / d parameters for a point lens model.

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            dA_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each epoch.
        """
        if self._magnification_objects is None:
            self._set_magnification_objects()

        d_A_d_params = {key: np.zeros(len(self.times)) for key in parameters}
        for method, selection in self.methods_indices.items():
            d_A_d_params_selection = \
                self._magnification_objects[method].get_d_A_d_params(parameters)
            for key in parameters:
                d_A_d_params[key][selection] = d_A_d_params_selection[key]

        return d_A_d_params

    def get_d_A_d_rho(self):
        """
        Calculate d A / d rho for a point lens model.

        No Inputs

        Returns :
            dA_drho: *np.array*
                Values are the partial derivatives for rho
                evaluated at each data point.
        """
        if self._magnification_objects is None:
            self._set_magnification_objects()

        d_A_d_rho = np.zeros(len(self.times))
        for method, selection in self.methods_indices.items():
            if method.lower() == 'point_source':
                d_A_d_rho_selection = np.zeros(np.sum(selection))
            else:
                d_A_d_rho_selection = \
                    self._magnification_objects[method].get_d_A_d_rho()

            d_A_d_rho[selection] = d_A_d_rho_selection

        return d_A_d_rho

    @property
    def methods_for_epochs(self):
        """
        *list*

        for each epoch, decide which methods should be used to
        calculate magnification, but don't run the calculations
        """
        if self._methods_for_epochs is None:
            out = [self._default_method] * len(self.times)
            if self._methods_epochs is None:
                return out

            brackets = np.searchsorted(self._methods_epochs, self.times)
            n_max = len(self._methods_epochs)

            out = [self._methods_names[value - 1]
                   if (value > 0 and value < n_max) else self._default_method
                   for value in brackets]
            self._methods_for_epochs = out

        return self._methods_for_epochs

    @property
    def methods_indices(self):
        """
        *dict*

        Keys are the magnification methods. Values are a boolean index that
        indicate which epochs should be calculated with each method.
        """
        if self._methods_indices is None:
            self._methods_indices = {}
            methods = self.methods_for_epochs
            methods_ = np.array(methods)

            for method in set(methods):
                selection = (methods_ == method)
                self._methods_indices[method] = selection

        return self._methods_indices
