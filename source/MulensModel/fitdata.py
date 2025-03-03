import numpy as np
import warnings

import MulensModel as mm


class FitData(object):
    """
    Performs a least squares linear fit for given dataset and model to
    determine the source flux(es) and (optionally) blend flux. After creating
    the object, you must run :py:func:`~update()` to perform the linear fit for
    the fluxes and calculate the chi2. To perform the linear fit without
    calculating chi2, you can run :py:func:`fit_fluxes()`. If you change
    anything in the object, e.g. the model parameters, you *must* re-run
    :py:func:`~update()` or :py:func:`~fit_fluxes()`.

    Arguments :
        model: :py:class:`~MulensModel.model.Model` object
            The model to fit to the data.

        dataset: :py:class:`~MulensModel.mulensdata.MulensData` object
            A single photometric dataset to be fitted.

        fix_blend_flux: *False* or *float*, optional
            Default is *False*, i.e. allow the blend flux to be a free
            parameter. If set to a float, it will fix the blend value to that
            value.

        fix_source_flux: *False*, *float*, or *list*, optional
            Default is *False*, i.e. allow the source flux to be a free
            parameter. If set to a float, it will fix the source value to that
            value. For binary source models, a list should be used to set the
            fluxes of the individual sources or fix one and not the other, e.g.
            [2.3, False] would fix source_flux_0 to 2.3 but allow a free fit to
            source_flux_1.

        fix_source_flux_ratio: *False* or *float*, optional
            For binary source models, source_flux_ratio is the flux ratio
            between two  components, i.e.,
            source_flux_ratio = source_flux_1 / source_flux_0
            Default is *False*, i.e. allow the source flux to be a free
            parameter. If set to a float, it will fix the source value to that
            value.

    """

    def __init__(self, model=None, dataset=None, fix_blend_flux=False,
                 fix_source_flux=False, fix_source_flux_ratio=False):
        self.model = model
        self.dataset = dataset

        # Setup limb-darkening
        self._gamma = 0.
        if self.model.parameters.is_finite_source():
            if self.dataset.bandpass is not None:
                try:
                    self._gamma = self.model.get_limb_coeff_gamma(
                        self.dataset.bandpass)
                except KeyError:
                    msg = (
                        'Dataset bandpass is {0} but model does not have a ' +
                        'limb-darkening coefficient for {0}. Assuming zero.')
                    warnings.warn(msg.format(self.dataset.bandpass))

        # fit parameters
        self.fix_blend_flux = fix_blend_flux
        self.fix_source_flux_ratio = self._set_fix_source_flux_ratio(fix_source_flux_ratio)
        self.fix_source_flux = self._set_fix_source_flux(fix_source_flux)

        # parameters fluxes of various sources
        self._source_fluxes = None
        self._blend_flux = None
        self._source_flux_ratio = None

        # chi2 parameters
        self._chi2_per_point = None
        self._chi2 = None
        self._chi2_gradient = None

        self._data_magnification = None
        self._data_magnification_curve = None
        self._data_magnification_curves = None

    def __getattr__(self, item):
        return object.__getattribute__(self, item)

    def _set_fix_source_flux(self, fix_source_flux):
        if fix_source_flux is False:
            return fix_source_flux
        else:
            if isinstance(fix_source_flux, list):
                if len(fix_source_flux) == self._model.n_sources:
                    return fix_source_flux

            elif isinstance(fix_source_flux, (float, int)) and (self._model.n_sources == 1):
                return [fix_source_flux]

        msg = ("you have {0}".format(self._model.n_sources) +
               " sources. Thus, fix_source_flux should be a list of" +
               "length {0}".format(self._model.n_sources) +
               "(or False).")
        raise ValueError(msg)

    def _set_fix_source_flux_ratio(self, fix_source_flux_ratio):
        if fix_source_flux_ratio is False:
            return fix_source_flux_ratio
        else:
            if isinstance(fix_source_flux_ratio, list):
                if len(fix_source_flux_ratio) == (self._model.n_sources - 1):
                    return fix_source_flux_ratio

            elif isinstance(fix_source_flux_ratio, (float, int)) and (self._model.n_sources == 2):
                return np.array([fix_source_flux_ratio])

        msg = ("you have {0}".format(self._model.n_sources) +
               " sources. Thus, fix_source_flux_ratio should be a list of" +
               "length {0} (or Flase).".format(self._model.n_sources - 1))
        raise ValueError(msg)

    def _check_for_flux_ratio_errors(self):
        """
        If combination of settings and models are invalid, raise exceptions.
        """

        if self.fix_source_flux_ratio is not False:
            if self._model.n_sources == 1:
                raise ValueError('fix_source_flux_ratio is not defined for only 1 source!')
            elif len(self.fix_source_flux_ratio) != self._model.n_sources-1:
                raise ValueError('fix_source_flux_ratio should be a *list* of len={0}'.format(self._model.n_sources-1))
            elif self.fix_source_flux is not False:
                msg = ('fix_source_flux_ratio + fixed_source_flux not ' +
                       'implemented. Fix the fluxes for each source ' +
                       'individually instead.')
                raise NotImplementedError(msg)

    def update(self, bad=False):
        """
        Calculate the best-fit source and blend fluxes as well as the chi2.

        Parameters :
            bad: *bool*
                Default is *False*. If *True* recalculates the data
                magnification for each point to ensure that there are values
                even for bad datapoints.

        No returns.
        """
        self.fit_fluxes()

        # Calculate chi2
        model_flux = self.get_model_fluxes(bad=bad)
        diff = self._dataset.flux - model_flux
        self._chi2_per_point = (diff / self._dataset.err_flux)**2

    def _set_data_magnification_curves(self, bad=True):
        if bad:
            select = np.ones(self._dataset.n_epochs, dtype=bool)
        else:
            select = self._dataset.good

        if self.dataset.ephemerides_file is None:
            satellite_skycoord = None
        else:
            satellite_skycoord = self.dataset.satellite_skycoord[select]

        magnification_kwargs = {
            'gamma': self.gamma, 'satellite_skycoord': satellite_skycoord}

        if self._model.n_sources == 1:
            self._data_magnification_curve = \
                self._model.get_magnification_curve(time=self._dataset.time[select], **magnification_kwargs)
        elif self._model.n_sources >= 2:
            self._data_magnification_curves = self._model.get_magnification_curves(
                        time=self._dataset.time[select], **magnification_kwargs)
            for i in range(self._model.n_sources):
                self.__setattr__('_data_magnification_curve_{0}'.format(i+1), self._data_magnification_curves[i])

    def _calculate_magnifications(self, bad=True):
        """
        Calculate the model magnifications for the epochs of the dataset.
        """
        self._set_data_magnification_curves(bad=bad)

        if self._model.n_sources == 1:
            mag_matrix = self._data_magnification_curve.get_magnification()
        elif self._model.n_sources >= 2:
            mag_matrix = []
            for i in range(self._model.n_sources):
                mag_matrix.append(self.__getattr__('_data_magnification_curve_{0}'.format(i+1)).get_magnification())

        else:
            msg = ("{0} ".format(self._model.n_sources) +
                   "sources used. Function model.get_magnification can only handle <=2 sources")
            raise NotImplementedError(msg)

        if bad:
            self._data_magnification = mag_matrix
        else:
            if self._model.n_sources == 1:
                self._data_magnification = np.zeros(self._dataset.n_epochs)
                self._data_magnification[self._dataset.good] = mag_matrix
            else:
                self._data_magnification = [np.zeros(self._dataset.n_epochs)]
                self._data_magnification[0][self._dataset.good] = mag_matrix[0]
                for source in range(1, self.model.n_sources):
                    self._data_magnification.append(np.zeros(self._dataset.n_epochs))
                    self._data_magnification[source][self._dataset.good] = mag_matrix[source]

    def _get_xy_qflux(self):
        """
        Apply a fixed flux ratio.
        flux = sum_i(f_i * A_i) + f_b
             = f_1 * [ A_1 + sum_i>1(q_i * A_i)] + f_b
        """
        y = self._dataset.flux[self._dataset.good]
        x = np.array(self._data_magnification[0][self._dataset.good])
        self.n_fluxes = 1
        for i in range(1, self._model.n_sources):
            if self.fix_source_flux_ratio[i-1] is False:
                x = np.vstack((x, self._data_magnification[i][self._dataset.good]))
                self.n_fluxes += 1
            else:
                if len(x.shape) == 1:
                    x += self.fix_source_flux_ratio[i-1] * self._data_magnification[i][self._dataset.good]
                else:
                    x[0, :] += self.fix_source_flux_ratio[i-1] * self._data_magnification[i][self._dataset.good]

        return (x, y)

    def _get_xy_individual_fluxes(self):
        """ Account for source fluxes individually """
        y = self._dataset.flux[self._dataset.good]

        if self.fix_source_flux is False:
            x = np.array(self._data_magnification)
            if self.model.n_sources == 1:
                x = x[self._dataset.good]
            else:
                x = x[:, self._dataset.good]

            self.n_fluxes = self._model.n_sources
        else:
            x = None
            if self._model.n_sources == 1:
                y -= self.fix_source_flux[0] * self._data_magnification[self._dataset.good]
            else:
                for i in range(self._model.n_sources):
                    if self.fix_source_flux[i] is False:
                        self.n_fluxes += 1
                        if x is None:
                            x = self._data_magnification[i][self._dataset.good]
                        else:
                            x = np.vstack((x, self._data_magnification[i][self._dataset.good]))

                    else:
                        y -= self.fix_source_flux[i] * self._data_magnification[i][self._dataset.good]

        return (x, y)

    def _setup_linalg_arrays(self):
        """
        Create xT and y arrays
        """
        (x, y) = self._create_arrays()
        xT = self._invert_x_array(x)
        (xT, y) = self._weight_linalg_arrays(xT, y)
        return (xT, y)

    def _create_arrays(self):
        """ Create x and y arrays"""
        # Initializations
        self.n_fluxes = 0
        n_epochs = np.sum(self._dataset.good)
        self._calculate_magnifications(bad=False)

        # Account for source fluxes
        if self.fix_source_flux_ratio is not False:
            self._check_for_flux_ratio_errors()
            (x, y) = self._get_xy_qflux()
        else:
            (x, y) = self._get_xy_individual_fluxes()

        # Account for free or fixed blending
        # Should do a runtime test to compare with lines 83-94
        if self.fix_blend_flux is False:
            self.n_fluxes += 1
            if x is None:
                x = np.ones((1, n_epochs))
            else:
                x = np.vstack((x, np.ones(n_epochs)))

        elif self.fix_blend_flux == 0.:
            pass
        else:
            y -= self.fix_blend_flux

        return (x, y)

    def _invert_x_array(self, x):
        """ Take the transpose of x """
        n_epochs = np.sum(self._dataset.good)
        xT = np.copy(x).T
        xT.shape = (n_epochs, self.n_fluxes)

        return xT

    def _weight_linalg_arrays(self, xT, y):
        """weight by data uncertainties"""
        # Take into account uncertainties
        sigma_inverse = 1. / self._dataset.err_flux[self._dataset.good]
        y *= sigma_inverse
        xT *= np.array([sigma_inverse] * self.n_fluxes).T

        return (xT, y)

    def fit_fluxes(self):
        """
        Execute the linear least squares fit to determine the fitted fluxes.
        Sets the values of :py:obj:`~source_fluxes`, :py:obj:`~blend_flux`,
        and (if applicable) :py:obj:`~source_flux`.

        Does *not* calculate chi2. To fit for the fluxes and calculate chi2,
        run :py:func:`~update()`.

        No parameters.

        No returns.
        """

        # Bypass this code if all fluxes are fixed.
        if isinstance(self.fix_source_flux, (list, float)):
            if isinstance(self.fix_blend_flux, (float)):
                proceed = False
                if isinstance(self.fix_source_flux, (list)):
                    for item in self.fix_source_flux:
                        if isinstance(item, (float)):
                            pass
                        else:
                            if item is False:
                                proceed = True

                if not proceed:
                    self._calculate_magnifications(bad=False)
                    self._blend_flux = self.fix_blend_flux
                    self._source_fluxes = np.array(self.fix_source_flux)
                    return

        (xT, y) = self._setup_linalg_arrays()

        # Solve for the coefficients in y = fs * x + fb (point source)
        # These values are: F_s1, F_s2,..., F_b.
        try:
            results = np.linalg.lstsq(xT, y, rcond=-1)[0]
        except ValueError as e:
            message = (
                "{0}\nIf either of these numbers ({1}, {2}) is greater than "
                "zero, there is a NaN somewhere, probably in the data. The "
                "cause of this error may be the epochs with extreme "
                "brightness (e.g., 99.999 mag), which is sometimes used to "
                "mark bad data. Other possible reason is mistakenly using "
                "phot_fmt='flux' instead of 'mag'")
            args = (e, np.sum(np.isnan(xT)), np.sum(np.isnan(y)))
            raise ValueError(message.format(*args))

        # Record the results
        if self.fix_source_flux_ratio is False:
            if self.fix_source_flux is False:
                source_fluxes = results[0:self._model.n_sources]
            else:
                source_fluxes = []
                index = 0
                for i in range(self._model.n_sources):
                    if self.fix_source_flux[i] is False:
                        source_fluxes.append(results[index])
                        index += 1
                    else:
                        source_fluxes.append(self.fix_source_flux[i])

        else:
            source_fluxes = results[0]
            j = 1
            for i in range(1, self._model.n_sources):
                if self.fix_source_flux_ratio[i-1] is False:
                    source_fluxes = np.hstack((source_fluxes, results[j]))
                    j += 1
                else:
                    source_fluxes = np.hstack((source_fluxes, results[0]*self.fix_source_flux_ratio[i-1]))

        self._source_fluxes = np.array(source_fluxes)

        if self.fix_blend_flux is False:
            self._blend_flux = results[-1]
        else:
            self._blend_flux = self.fix_blend_flux

    def get_data_magnification(self, bad=False):
        """
        Calculates the model magnification for each data point.

        Parameters :
            bad: *boolean*
                If *True*, calculates the magnification for all points.
                If *False*, only calculates the magnification for good data
                points. Values for bad data points are set to 0. Default is
                *False*.

        Returns :
            data_magnification: *np.ndarray*
                The model magnification evaluated for each datapoint. If there
                is more than one source, the magnification of each source is
                reported separately.
        """

        self._calculate_magnifications(bad=bad)
        return self._data_magnification

    def get_model_fluxes(self, bad=False):
        """
        Calculate model in flux space.

        Parameters :
            bad: *bool*
                Default is *False*. If *True* recalculates the data
                magnification for each point to ensure that the values
                for bad datapoints are calculated (otherwise, they are set to
                the magnitude of the blend).

        Returns :
            model_flux: *np.ndarray*
                The model flux evaluated for each datapoint.
        """
        if self.source_fluxes is None:
            raise AttributeError('you need to run FitData.fit_fluxes() first to execute the linear fit.')

        if bad:
            self._calculate_magnifications(bad=True)

        return self._get_model_flux(self.source_fluxes, self.blend_flux)

    def _get_model_flux(self, source_fluxes, blend_flux):
        """
        Calculate model flux based on self._data_magnification.
        """
        model_flux = blend_flux * np.ones(self._dataset.n_epochs)

        if self._model.n_sources == 1:
            model_flux += source_fluxes * self._data_magnification
        else:
            for i in range(self._model.n_sources):
                model_flux += source_fluxes[i] * self._data_magnification[i]

        return model_flux

    def get_model_magnitudes(self, **kwargs):
        """
        Calculate model in magnitude space

        Parameters :
            ``**kwargs``:
                see :py:func:`get_model_fluxes()`

        Returns :
            model_mag: *np.ndarray*
                The model magnitude evaluated for each datapoint.
        """
        model_flux = self.get_model_fluxes(**kwargs)
        model_mag = mm.Utils.get_mag_from_flux(model_flux)

        return model_mag

    def scale_fluxes(self, source_flux, blend_flux):
        """
        Rescale the data fluxes to an arbitrary flux scale:
            flux = source_flux_0 * (data.flux - blend_flux) / source_flux
            flux += blend_flux_0
            err_flux = source_flux_0 * data.err_flux / source_flux

        Parameters :
            source_flux: *float*, *list*, *np.array*
                Flux of the source in the desired system. If n_sources > 1 and
                source_flux has more than one element, the elements are
                summed to produce the overall scaling flux.

            blend_flux: *float*
                Flux of the blend in the desired system

        Returns :
            flux: *np.ndarray*
                Fluxes from the data rescaled to the desired system.

            err_flux: *np.ndarray*
                Uncertainties of fluxes from the data rescaled to the desired
                system.
        """
        if self.model.n_sources == 1:
            data_source_flux = self.source_flux
        else:
            data_source_flux = np.sum(self.source_fluxes)
            if len(source_flux) > 1:
                source_flux = np.sum(source_flux)

        flux = source_flux * (self._dataset.flux - self.blend_flux)
        flux /= data_source_flux
        flux += blend_flux
        err_flux = source_flux * self._dataset.err_flux / data_source_flux

        return (flux, err_flux)

    def get_residuals(self, phot_fmt=None, source_flux=None, blend_flux=None,
                      bad=False):
        """
        Calculate the residuals for each datapoint relative to the model.

        Parameters :
            phot_fmt: *str*, optional
                specify whether the residuals should be returned in
                magnitudes ('mag') or in flux ('flux'). Default is
                'mag'. If 'scaled', will return the residuals in magnitudes
                scaled to source_flux, blend_flux.

            source_flux, blend_flux: *float*
                reference source and blend fluxes for scaling the residuals

            bad: *bool*
                Default is *False*. If *True* recalculates the data
                magnification for each point to ensure that there are values
                even for bad datapoints.

        Returns :
            residuals: *np.ndarray*
                the residuals for the corresponding dataset.

            errorbars: *np.ndarray*
                the scaled errorbars for each point. For plotting
                errorbars for the residuals.
        """
        if bad:
            self._calculate_magnifications(bad=True)

        if phot_fmt == 'mag':
            warnings.warn(
                '"mag" returns residuals in the original data flux system. To scale the residuals, use "scaled".')
            residuals = self._dataset.mag - self.get_model_magnitudes()
            errorbars = self._dataset.err_mag
        elif phot_fmt == 'flux':
            residuals = self._dataset.flux - self.get_model_fluxes()
            errorbars = self._dataset.err_flux
        elif phot_fmt == 'scaled':
            if source_flux is None or blend_flux is None:
                raise ValueError('If phot_fmt=scaled, source_flux and blend_flux must also be specified.')

            model_flux = self._get_model_flux(source_flux, blend_flux)
            model_mag = mm.Utils.get_mag_from_flux(model_flux)
            (flux, err_flux) = self.scale_fluxes(source_flux, blend_flux)
            (mag, errorbars) = mm.Utils.get_mag_and_err_from_flux(flux, err_flux)
            residuals = mag - model_mag
        else:
            raise ValueError(
                'phot_fmt must be one of "mag", "flux", or "scaled". Your value: {0}'.format(phot_fmt))

        return (residuals, errorbars)

    def _check_for_gradient_implementation(self, parameters):
        """
        Check that the gradient methods are implemented for the requested
        values.
        """
        # Implemented for the requested parameters?
        if not isinstance(parameters, list):
            parameters = [parameters]
        implemented = {'t_0', 't_E', 'u_0', 't_eff', 'pi_E_N', 'pi_E_E', 'rho'}
        if len(set(parameters) - implemented) > 0:
            raise NotImplementedError((
                "chi^2 gradient is implemented only for {:}\nCannot work with {:}").format(implemented, parameters))

        # Implemented for the number of lenses in the model?
        if self.model.n_lenses != 1:
            raise NotImplementedError('chi2_gradient() only implemented for single lens models')

        if self.model.parameters.is_xallarap:
            raise NotImplementedError('Gradient for xallarap models is not implemented yet')

    def get_chi2_gradient(self, parameters):
        """
        Fits fluxes and calculates chi^2 gradient (also called Jacobian), i.e.,
        :math:`d chi^2/d parameter`.

        Parameters :
            parameters: *str* or *list*, required
                Parameters with respect to which gradient is calculated.
                Currently accepted parameters are: ``t_0``, ``u_0``, ``t_eff``,
                ``t_E``, ``pi_E_N``, and ``pi_E_E``. The parameters for
                which you request gradient must be defined in py:attr:`~model`.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient
        """
        self.fit_fluxes()
        self.calculate_chi2_gradient(parameters)
        return self.chi2_gradient

    def calculate_chi2_gradient(self, parameters):
        """
        Calculates chi^2 gradient (also called Jacobian), i.e.,
        :math:`d chi^2/d parameter` WITHOUT refitting for the fluxes. Saves
        computations if, e.g., you want to retrieve both py:attr:`~chi2` and
        py:attr:`~chi2_gradient`.

        Parameters :
            parameters: *str* or *list*, required
                Parameters with respect to which gradient is calculated.
                Currently accepted parameters are: ``t_0``, ``u_0``, ``t_eff``,
                ``t_E``, ``pi_E_N``, and ``pi_E_E``. The parameters for
                which you request gradient must be defined in py:attr:`~model`.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient
        """
        self._check_for_gradient_implementation(parameters)

        # Calculate factor
        flux_factor = self.get_model_fluxes() - self.dataset.flux
        flux_factor *= 2. * self.source_flux / self.dataset.err_flux**2

        gradient = self.get_d_A_d_params_for_point_lens_model(parameters)
        for (key, value) in gradient.items():
            gradient[key] = np.sum((flux_factor[self.dataset.good] * value))

        if len(parameters) == 1:
            out = gradient[parameters[0]]
        else:
            out = np.array([gradient[p] for p in parameters])

        self._chi2_gradient = out

        return self._chi2_gradient

    def get_d_A_d_params_for_point_lens_model(self, parameters):
        """
        Calculate d A / d parameters for a point lens model.

        Parameters :
            parameters: *list*
                List of the parameters to take derivatives with respect to.

        Returns :
            dA_dparam: *dict*
                Keys are parameter names from *parameters* argument above.
                Values are the partial derivatives for that parameter
                evaluated at each data point.
        """
        if self._model.n_sources > 1:
            raise NotImplementedError("gradient for multi-source models is not implemented/tested")

        if self._data_magnification_curve is None:
            self._set_data_magnification_curves()

        d_A_d_params = self._data_magnification_curve.get_d_A_d_params(parameters)

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
        if self._model.n_sources > 1:
            raise NotImplementedError("gradient for multi-source models is not implemented/tested")

        if 'rho' not in self.model.parameters.parameters:
            raise AttributeError('dA/drho cannot be calculated for a model without rho')

        if self._data_magnification_curve is None:
            self._set_data_magnification_curves()

        d_A_d_params = self._data_magnification_curve.get_d_A_d_rho()

        return d_A_d_params

    def get_dataset_trajectory(self):
        """
        Retrieve a :py:class:`~MulensModel.trajectory.Trajectory` object. If
        the :py:attr:`~dataset` has an ephemerides_file, apply it to the
        Trajectory, even if it is not part of the :py:attr:`~model`.

        No parameters.

        Returns :
            trajectory: :py:class:`~MulensModel.trajectory.Trajectory`
                Trajectory for given dataset.
        """
        kwargs = dict()
        if self.dataset.ephemerides_file is not None:
            kwargs['satellite_skycoord'] = self.dataset.satellite_skycoord

        return self.model.get_trajectory(self.dataset.time, **kwargs)

    def get_d_A_d_u_for_PSPL_model(self):
        raise NotImplementedError('This function was deprecated in Version 3.')

    def get_d_A_d_u_for_FSPL_model(self):
        raise NotImplementedError('This function was deprecated in Version 3.')

    def get_d_A_d_u_for_point_lens_model(self):
        raise NotImplementedError('This function was deprecated in Version 3.')

    @property
    def chi2_gradient(self):
        """
        *float* or *np.ndarray*

        Previously calculated chi^2 gradient (also called Jacobian),
        i.e., :math:`d chi^2/d parameter`. See :py:func:`~get_chi2_gradient()`
        and :py:func:`~calculate_chi2_gradient()`.

        Gives *None* if the chi2 gradient was not
        previously calculated using one of the functions mentioned
        above.
        """
        try:
            return self._chi2_gradient
        except AttributeError:
            return None

    @property
    def chi2(self):
        """
        *float*
        The total chi2 for the fitted dataset. Good points only.
        See :py:obj:`~MulensModel.mulensdata.MulensData.good`.

        If *None*, you need to run :py:func:`~update()` to execute the
        linear fit and calculate the chi2.
        """
        if self.chi2_per_point is None:
            return None
        else:
            return np.sum(self.chi2_per_point[self._dataset.good])

    @property
    def chi2_per_point(self):
        """
        *np.ndarray*

        The chi^2 contribution from each data point,
        e.g., ``chi2_per_point[k]`` returns the chi2 contribution
        from the *k*-th point of :py:obj:`dataset`. Includes bad
        datapoints.

        If *None*, you need to run :py:func:`~update()` to execute
        the linear fit and calculate the chi2.
        """
        return self._chi2_per_point

    @property
    def source_flux(self):
        """
        *float*

        The fitted source flux. Only defined for models with a single
        source. See also :py:obj:`~source_fluxes`

        If *None*, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        if self._model.n_sources == 1:
            return self.source_fluxes[0]
        else:
            msg = ("source_flux is defined only for models" +
                   " with ONE source, you have" +
                   " {0}".format(self._model.n_sources) +
                   " sources. Try FitData.source_fluxes instead")

            raise NameError(msg)

    @property
    def source_fluxes(self):
        """
        *np.array*

        The fitted source flux(es).

        If *None*, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        return self._source_fluxes

    @property
    def blend_flux(self):
        """
        *float*

        The fitted blend flux or the value set by
        fix_blend_flux (see :ref:`keywords`).

        If *None*, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        return self._blend_flux

    @property
    def source_flux_ratio(self):
        """
        *float*

        source_flux_ratio = source_flux_1 / source_flux_0

        i.e., the ratio of the fitted source fluxes or the value set by
        fix_source_flux_ratio (see :ref:`keywords`).

        If *None*, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        if self._model.n_sources == 1:
            msg = ("source_flux is defined only for models" +
                   " with multiple sources; you have 1 source.")
            raise NameError(msg)
        else:
            source_flux_ratios = []
            for i in range(1, self._model.n_sources):
                if (self.fix_source_flux_ratio is False) or (self.fix_source_flux_ratio[i-1] is False):
                    source_flux_ratios.append(self.source_fluxes[i] / self.source_fluxes[0])
                else:
                    source_flux_ratios.append(self.fix_source_flux_ratio[i-1])

            return source_flux_ratios

    @property
    def dataset(self):
        """
        :py:class:`~MulensModel.mulensdata.MulensData`

        A single photometric dataset to be fitted.
        """
        return self._dataset

    @dataset.setter
    def dataset(self, new_value):
        if not isinstance(new_value, mm.MulensData):
            raise TypeError("Dataset has to of MulensData type, not: " +
                            str(type(new_value)))
        self._dataset = new_value

    @property
    def model(self):
        """
        :py:class:`~MulensModel.model.Model`

        The model to fit to the data.
        """
        return self._model

    @model.setter
    def model(self, new_value):
        self._model = new_value

    @property
    def data_magnification(self):
        """
        Returns previously calculated magnifications. To calculate the
        magnifications (e.g., if something changed in the model), use
        :py:func:`~get_data_magnification()`.

        Returns :
            data_magnification: *np.ndarray*
                The model magnification evaluated for each datapoint. If there
                is more than one source, the magnification of each source is
                reported separately.
        """
        return self._data_magnification

    @property
    def magnification_curve(self):
        """
        Returns previously calculated magnification curve.
        """
        return self._data_magnification_curve

    @property
    def magnification_curves(self):
        """
        Returns previously calculated magnification curves.

        Returns :
            *tuple* of *:py:class:`~MulensModel.magnification.MagnificationCurve* objects
            i.e., the model magnification curve evaluated for each datapoint.
        """
        out = ()
        for i in range(self._model.n_sources):
            out.append(self.__getattr__('_data_magnification_curve_{0}'.format(i+1)).get_magnification())

        return out

    @property
    def gamma(self):
        """
        *float*

        Limb-darkening coefficient for this fit. Set by
        :py:attr:`~MulensModel.mulensdata.MulensData.bandpass` and
        :py:func:`~MulensModel.model.Model.get_limb_coeff_gamma()`.
        """
        return self._gamma

    class FSPL_Derivatives(object):

        def __init__(self, fit):
            raise NotImplementedError(
                'The FSPL_Derivatives class was deprecated in Version 3. ' +
                'Its various functions were incorporated into the new PointLens classes.')
