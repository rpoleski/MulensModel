import numpy as np
import warnings
from MulensModel.trajectory import Trajectory
import astropy.units as u


class FitData:
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
            [2.3, False] would fix f_source_0 to 2.3 but allow a free fit to
            f_source_1.

        fix_q_flux: *False* or *float*, optional
            For binary source models, q_flux is the flux ratio between two
            components, i.e. q_flux = f_source_1 / f_source_0
            Default is *False*, i.e. allow the source flux to be a free
            parameter. If set to a float, it will fix the source value to that
            value.

    """

    def __init__(self, model=None, dataset=None, fix_blend_flux=False,
                 fix_source_flux=False, fix_q_flux=False):
        self.model = model
        self.dataset = dataset

        # Setup limb-darkening
        if self.dataset.bandpass is None:
            self.gamma = 0.
        else:
            try:
                self.gamma = self.model.get_limb_coeff_gamma(
                    self.dataset.bandpass)
            except KeyError:
                warnings.warn(
                    'Dataset bandpass is {0} '.format(self.dataset.bandpass) +
                    'but model does not have a limb-darkening coefficient ' +
                    'for {0}. Assuming zero.'.format(self.dataset.bandpass))
                self.gamma = 0.

        # fit parameters
        self.fix_blend_flux = fix_blend_flux
        self.fix_q_flux = fix_q_flux
        if isinstance(fix_source_flux, list) or (fix_source_flux is False):
            self.fix_source_flux = fix_source_flux
        else:
            if self._model.n_sources == 1:
                self.fix_source_flux = [fix_source_flux]
            else:
                msg = ("you have {0}".format(self._model.n_sources) +
                       " sources. Thus, fix_source_flux should be a list of" +
                       "length {0}".format(self._model.n_sources) +
                       "(or False).")
                raise ValueError(msg)

        # parameters fluxes of various sources
        self._source_fluxes = None
        self._blend_flux = None
        self._q_flux = None

        # chi2 parameters
        self._chi2_per_point = None
        self._chi2 = None

    def _check_for_q_flux_errors(self):
        """
        If combination of settings and models are invalid, raise exceptions.
        """

        if self.fix_q_flux is not False:
            if self._model.n_sources != 2:
                msg = ('fix_q_flux only valid for models with 2 sources.' +
                       'n_sources = {0}'.format(self._model.n_sources))
                raise ValueError(msg)
            elif self.fix_source_flux is not False:
                msg = ('fix_q_flux + fixed_source_flux not implemented.' +
                       'Fix the fluxes for each source individually instead.')
                raise NotImplementedError(msg)

    def update(self):
        """
        Calculate the best-fit source and blend fluxes as well as the chi2.
        """
        self.fit_fluxes()

        # Calculate chi2
        model_flux = self.get_model_fluxes()
        diff = self._dataset.flux - model_flux
        self._chi2_per_point = (diff / self._dataset.err_flux)**2

    def _calc_magnifications(self):
        """
        Calculate the model magnifications for the good epochs of the dataset.
        """
        # select = self._dataset.good # Is this faster?

        # currently, model.magnification is good for up to two
        # sources
        if self._model.n_sources == 1:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[self._dataset.good], gamma=self.gamma)
        elif self._model.n_sources == 2:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[self._dataset.good], gamma=self.gamma,
                separate=True)
        else:
            msg = ("{0}".format(self._model.n_sources) +
                   " sources used. model.magnification can only" +
                   " handle <=2 sources")
            raise NotImplementedError(msg)

        self._data_magnification = mag_matrix

    def _get_xy_qflux(self):
        """Use a flux ratio constraint"""
        y = self._dataset.flux[self._dataset.good]
        x = np.array(
            self._data_magnification[0] +
            self.fix_q_flux * self._data_magnification[1])
        self.n_fluxes = 1

        return (x, y)

    def _get_xy_individual_fluxes(self):
        """ Account for source fluxes individually """
        y = self._dataset.flux[self._dataset.good]

        if self.fix_source_flux is False:
            x = np.array(self._data_magnification)
            self.n_fluxes = self._model.n_sources
        else:
            x = None
            if self._model.n_sources == 1:
                y -= self.fix_source_flux[0] * self._data_magnification
            else:
                for i in range(self._model.n_sources):
                    if self.fix_source_flux[i] is False:
                        self.n_fluxes += 1
                        if x is None:
                            if self._model.n_sources == 1:
                                x = np.array(self._data_magnification)
                            else:
                                x = np.array(self._data_magnification[i])

                        else:
                            x = np.vstack((x, self._data_magnification[i]))

                    else:
                        y -= (self.fix_source_flux[i] *
                              self._data_magnification[i])

        return (x, y)

    def _setup_linalg_arrays(self):
        """
        Create xT and y arrays
        """
        # Initializations
        self.n_fluxes = 0
        n_epochs = np.sum(self._dataset.good)
        self._calc_magnifications()

        # Account for source fluxes
        if self.fix_q_flux is not False:
            self._check_for_q_flux_errors()
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

        # Take the transpose of x and weight by data uncertainties
        xT = np.copy(x).T
        xT.shape = (n_epochs, self.n_fluxes)

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
        """

        (xT, y) = self._setup_linalg_arrays()

        # Solve for the coefficients in y = fs * x + fb (point source)
        # These values are: F_s1, F_s2,..., F_b.
        try:
            results = np.linalg.lstsq(xT, y, rcond=-1)[0]
        except ValueError as e:
            raise ValueError(
                '{0}\n'.format(e) +
                'If either of these numbers ({0}, {1})'.format(
                    np.sum(np.isnan(xT)), np.sum(np.isnan(y))) +
                ' is greater than zero, there is a NaN somewhere,' +
                ' probably in the data.')

        # Record the results
        if self.fix_q_flux is False:
            if self.fix_source_flux is False:
                self._source_fluxes = results[0:self._model.n_sources]
            else:
                self._source_fluxes = []
                index = 0
                for i in range(self._model.n_sources):
                    if self.fix_source_flux[i] is False:
                        self._source_fluxes.append(results[index])
                        index += 1
                    else:
                        self._source_fluxes.append(self.fix_source_flux[i])

        else:
            self._source_fluxes = [results[0], results[0] * self.fix_q_flux]

        if self.fix_blend_flux is False:
            self._blend_flux = results[-1]
        else:
            self._blend_flux = self.fix_blend_flux

    def get_model_fluxes(self):
        """
        Returns :
            model_flux: *np.ndarray*
                The model flux evaluated for each datapoint.
        """
        if self.source_fluxes is None:
            raise AttributeError(
                'you need to run FitData.fit_fluxes() first to execute the' +
                'linear fit.')

        model_flux = np.ones(self._dataset.n_epochs) * self.blend_flux
        if self._model.n_sources == 1:
            model_flux += self.source_flux * self._data_magnification
        else:
            for i in range(self._model.n_sources):
                model_flux += self.source_fluxes[i] \
                        * self._data_magnification[i]

        return model_flux

    def _check_for_gradient_implementation(self, parameters):
        """ Check that the gradient methods are implemented for the requested
        values. """

        # Implemented for the requested parameters?
        if not isinstance(parameters, list):
            parameters = [parameters]
        implemented = {'t_0', 't_E', 'u_0', 't_eff', 'pi_E_N', 'pi_E_E'}
        if len(set(parameters) - implemented) > 0:
            raise NotImplementedError((
                "chi^2 gradient is implemented only for {:}\nCannot work " +
                "with {:}").format(implemented, parameters))
        gradient = {param: 0 for param in parameters}

        # Implemented for the number of sources in the model?
        if self.model.n_lenses != 1:
            raise NotImplementedError(
                'chi2_gradient() works only single lens models currently')

        # Implemented for finite source effects?
        if 'rho' in parameters or 't_star' in parameters:
            as_dict = self.model.parameters.as_dict()
            if 'rho' in as_dict or 't_star' in as_dict:
                raise NotImplementedError(
                    'Event.chi2_gradient() is not working ' +
                    'for finite source models yet')

    def get_chi2_gradient(self, parameters):
        """ Same as :py:func:`~chi2_gradient`, but fits for the fluxes first."""
        self.fit_fluxes()
        return self.chi2_gradient()

    def chi2_gradient(self, parameters):
        """
        Calculate chi^2 gradient (also called Jacobian), i.e.,
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
        self._check_for_gradient_implementation(parameters)

        # Setup
        gradient = {param: 0 for param in parameters}
        as_dict = self.model.parameters.as_dict()

        # ***seems unnecessary
        # # Define a Fit given the model and perform linear fit for fs and fb
        # self._update_data_in_model()
        # self.fit = Fit(
        #     data=self.datasets, magnification=self.model.data_magnification)
        # # For binary source cases, the above line would need to be replaced,
        # # so that it uses self.model.fit.
        # if fit_blending is not None:
        #     self.fit.fit_fluxes(fit_blending=fit_blending)
        # else:
        #     self.fit.fit_fluxes()
        #

        # ***
        # JCY - Originally implemented for arbitrary chi2 fmt, but this should
        # always be fluxes, right?

        # for (i, dataset) in enumerate(self.datasets):
        #     (data, err_data) = dataset.data_and_err_in_chi2_fmt()
        #     factor = data - self.fit.get_chi2_format(data=dataset)
        #     factor *= -2. / err_data**2
        #     if dataset.chi2_fmt == 'mag':
        #         factor *= -2.5 / (log(10.) * Utils.get_flux_from_mag(data))
        #     factor *= self.fit.flux_of_sources(dataset)[0]

        # Calculate factor
        # JCY - Everything below here should be refactored into smaller bits.
        factor = self.dataset.flux - self.get_model_fluxes()
        factor *= -2. / self.dataset.err_flux**2
        factor *= self.source_flux

        # Get source location
        trajectory = self.model.get_trajectory(self.dataset.time)
        u_2 = trajectory.x**2 + trajectory.y**2
        u_ = np.sqrt(u_2)

        # Calculate derivatives
        d_A_d_u = -8. / (u_2 * (u_2 + 4) * np.sqrt(u_2 + 4))
        factor *= d_A_d_u
        factor_d_x_d_u = (factor * trajectory.x / u_)[self.dataset.good]
        sum_d_x_d_u = np.sum(factor_d_x_d_u)
        factor_d_y_d_u = (factor * trajectory.y / u_)[self.dataset.good]
        sum_d_y_d_u = np.sum(factor_d_y_d_u)
        dt = self.dataset.time[self.dataset.good] - as_dict['t_0']

        # Exactly 2 out of (u_0, t_E, t_eff) must be defined and
        # gradient depends on which ones are defined.
        if 't_eff' not in as_dict:
            t_E = as_dict['t_E'].to(u.day).value
            if 't_0' in parameters:
                gradient['t_0'] += -sum_d_x_d_u / t_E
            if 'u_0' in parameters:
                gradient['u_0'] += sum_d_y_d_u
            if 't_E' in parameters:
                gradient['t_E'] += np.sum(factor_d_x_d_u * -dt / t_E**2)
        elif 't_E' not in as_dict:
            t_eff = as_dict['t_eff'].to(u.day).value
            if 't_0' in parameters:
                gradient['t_0'] += -sum_d_x_d_u * as_dict['u_0'] / t_eff
            if 'u_0' in parameters:
                gradient['u_0'] += sum_d_y_d_u + np.sum(
                        factor_d_x_d_u * dt / t_eff)
            if 't_eff' in parameters:
                gradient['t_eff'] += np.sum(
                        factor_d_x_d_u * -dt *
                        as_dict['u_0'] / t_eff**2)
        elif 'u_0' not in as_dict:
            t_E = as_dict['t_E'].to(u.day).value
            t_eff = as_dict['t_eff'].to(u.day).value
            if 't_0' in parameters:
                gradient['t_0'] += -sum_d_x_d_u / t_E
            if 't_E' in parameters:
                gradient['t_E'] += (
                        np.sum(factor_d_x_d_u * dt) -
                        sum_d_y_d_u * t_eff) / t_E**2
            if 't_eff' in parameters:
                gradient['t_eff'] += sum_d_y_d_u / t_E
        else:
            raise KeyError(
                'Something is wrong with ModelParameters in ' +
                'Event.chi2_gradient():\n', as_dict)

        # Below we deal with parallax only.
        if 'pi_E_N' in parameters or 'pi_E_E' in parameters:
            parallax = {
                'earth_orbital': False,
                'satellite': False,
                'topocentric': False}
            # JCY Not happy about this as it requires importing from other
            # modules. It is inelegant, which in my experience often means it
            # needs to be refactored.
            kwargs = {}
            if self.dataset.ephemerides_file is not None:
                kwargs['satellite_skycoord'] = self.dataset.satellite_skycoord
                
            trajectory_no_piE = Trajectory(
                self.dataset.time, self.model.parameters, parallax,
                self.model.coords, **kwargs)
            dx = (trajectory.x - trajectory_no_piE.x)[self.dataset.good]
            dy = (trajectory.y - trajectory_no_piE.y)[self.dataset.good]
            delta_E = dx * as_dict['pi_E_E'] + dy * as_dict['pi_E_N']
            delta_N = dx * as_dict['pi_E_N'] - dy * as_dict['pi_E_E']
            det = as_dict['pi_E_N']**2 + as_dict['pi_E_E']**2

            if 'pi_E_N' in parameters:
                gradient['pi_E_N'] += np.sum(
                    factor_d_x_d_u * delta_N + factor_d_y_d_u * delta_E)
                gradient['pi_E_N'] /= det
            if 'pi_E_E' in parameters:
                gradient['pi_E_E'] += np.sum(
                    factor_d_x_d_u * delta_E - factor_d_y_d_u * delta_N)
                gradient['pi_E_E'] /= det

        if len(parameters) == 1:
            out = gradient[parameters[0]]
        else:
            out = np.array([gradient[p] for p in parameters])
            
        return out

    @property
    def chi2(self):
        """
        Returns :
            chi2: *float*
                the total chi2 for the fitted dataset. Good points only. See
                :py:obj:`~MulensModel.mulensdata.MulensData.good`.

        If None, you need to run :py:func:`~update()` to execute the
        linear fit and calculate the chi2.
        """
        return np.sum(self.chi2_per_point[self._dataset.good])

    @property
    def chi2_per_point(self):
        """
        Returns :
            chi2_per_point: *np.ndarray*
                Chi^2 contribution from each data point,
                e.g. ``chi2_per_point[k]`` returns the chi2 contribution
                from the *k*-th point of :py:obj:`dataset`. Includes bad
                datapoints.

        If None, you need to run :py:func:`~update()` to execute the
        linear fit and calculate the chi2.
        """
        return self._chi2_per_point

    @property
    def source_flux(self):
        """
        Returns :
            source_flux: *float*
            the fitted source flux. Only defined for models with a single
            source. See also :py:obj:`~source_fluxes`

        If None, you need to run :py:func:`~fit_fluxes()` or
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
        Returns :
            source_fluxes: *np.array*
                the fitted source flux(es).

        If None, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        return self._source_fluxes

    @property
    def blend_flux(self):
        """
        Returns :
            blend_flux: *float*
                the fitted blend flux or the value set by
                fix_blend_flux (see :ref:`keywords`).

        If None, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        return self._blend_flux

    @property
    def q_flux(self):
        """
        q_flux = f_source_1 / f_source_0

        Returns :
            q_flux: *float*
                the ratio of the fitted source fluxes or the value set by
                fix_q_flux (see :ref:`keywords`).

        If None, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        if self._model.n_sources != 2:
            msg = ("source_flux is defined only for models" +
                   " with TWO sources, you have" +
                   " {0}".format(self._model.n_sources) +
                   " sources.")
            raise NameError(msg)

        if self.fix_q_flux:
            return self.fix_q_flux
        else:
            return self.source_fluxes[1] / self.source_fluxes[0]

    @property
    def dataset(self):
        """
        :py:class:`~MulensModel.mulensdata.MulensData` object

        A single photometric dataset to be fitted.
        """
        return self._dataset

    @dataset.setter
    def dataset(self, new_value):
        self._dataset = new_value

    @property
    def model(self):
        """
        :py:class:`~MulensModel.model.Model` object

        The model to fit to the data.
        """
        return self._model

    @model.setter
    def model(self, new_value):
        self._model = new_value
