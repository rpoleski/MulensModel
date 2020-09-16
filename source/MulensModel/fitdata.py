import numpy as np
import warnings
from MulensModel.trajectory import Trajectory
from MulensModel.satelliteskycoord import SatelliteSkyCoord
from MulensModel.utils import Utils
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
        self.fix_source_flux_ratio = fix_source_flux_ratio
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
        self._source_flux_ratio = None

        # chi2 parameters
        self._chi2_per_point = None
        self._chi2 = None

    def _check_for_flux_ratio_errors(self):
        """
        If combination of settings and models are invalid, raise exceptions.
        """

        if self.fix_source_flux_ratio is not False:
            if self._model.n_sources != 2:
                msg = ('fix_source_flux_ratio only valid for models with 2 sources.' +
                       'n_sources = {0}'.format(self._model.n_sources))
                raise ValueError(msg)
            elif self.fix_source_flux is not False:
                msg = ('fix_source_flux_ratio + fixed_source_flux not implemented.' +
                       'Fix the fluxes for each source individually instead.')
                raise NotImplementedError(msg)

    def update(self, bad=False):
        """
        Calculate the best-fit source and blend fluxes as well as the chi2.

        Keywords :
            bad: *bool*
                Default is *False*. If *True* recalculates the data
                magnification for each point to ensure that there are values
                even for bad datapoints.

        """
        self.fit_fluxes()

        # Calculate chi2
        model_flux = self.get_model_fluxes(bad=bad)
        diff = self._dataset.flux - model_flux
        self._chi2_per_point = (diff / self._dataset.err_flux)**2

    def _calc_magnifications(self, bad=True):
        """
        Calculate the model magnifications for the epochs of the dataset.
        """
        if bad:
            select = np.ones(self._dataset.n_epochs, dtype=bool)
        else:
            select = self._dataset.good

        if self.dataset.ephemerides_file is None:
            satellite_skycoord = None
        else:
            satellite_skycoords = SatelliteSkyCoord(
                 ephemerides_file=self.dataset.ephemerides_file)
            satellite_skycoord = satellite_skycoords.get_satellite_coords(
                self.dataset.time)

        magnification_kwargs = {
            'gamma': self.gamma, 'satellite_skycoord': satellite_skycoord}

        # currently, model.magnification is good for up to two
        # sources
        if self._model.n_sources == 1:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[select],
                **magnification_kwargs)
        elif self._model.n_sources == 2:
            mag_matrix = self._model.magnification(
                time=self._dataset.time[select], separate=True,
                **magnification_kwargs)
        else:
            msg = ("{0}".format(self._model.n_sources) +
                   " sources used. model.magnification can only" +
                   " handle <=2 sources")
            raise NotImplementedError(msg)

        if bad:
            self._data_magnification = mag_matrix
        else:
            if self._model.n_sources == 1:
                self._data_magnification = np.zeros(
                    self._dataset.n_epochs)
                self._data_magnification[self._dataset.good] = mag_matrix
            else:
                self._data_magnification = [np.zeros(self._dataset.n_epochs)]
                self._data_magnification[0][self._dataset.good] = mag_matrix[0]
                for source in range(1, self.model.n_sources):
                    self._data_magnification.append(
                        np.zeros(self._dataset.n_epochs))
                    self._data_magnification[
                        source][self._dataset.good] = mag_matrix[source]

    def _get_xy_qflux(self):
        """ Apply a fixed flux ratio. """
        y = self._dataset.flux[self._dataset.good]
        x = np.array(
            self._data_magnification[0][self._dataset.good] +
            self.fix_source_flux_ratio * self._data_magnification[1][self._dataset.good])
        self.n_fluxes = 1

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
                y -= (self.fix_source_flux[0] *
                    self._data_magnification[self._dataset.good])
            else:
                for i in range(self._model.n_sources):
                    if self.fix_source_flux[i] is False:
                        self.n_fluxes += 1
                        if x is None:
                            x = self._data_magnification[i][self._dataset.good]
                        else:
                            x = np.vstack(
                                (x, self._data_magnification[i][
                                    self._dataset.good]))

                    else:
                        y -= (self.fix_source_flux[i] *
                              self._data_magnification[i][self._dataset.good])

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
        # Initializations
        self.n_fluxes = 0
        n_epochs = np.sum(self._dataset.good)
        self._calc_magnifications(bad=False)

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
        if self.fix_source_flux_ratio is False:
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
            self._source_fluxes = [results[0], results[0] * self.fix_source_flux_ratio]

        if self.fix_blend_flux is False:
            self._blend_flux = results[-1]
        else:
            self._blend_flux = self.fix_blend_flux

    def get_data_magnification(self, bad=False):
        """
        Calculates the model magnification for each data point.

        Arguments :
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

        self._calc_magnifications(bad=bad)
        return self._data_magnification

    def get_model_fluxes(self, bad=False):
        """

        Keywords :
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
            raise AttributeError(
                'you need to run FitData.fit_fluxes() first to execute the' +
                'linear fit.')

        if bad:
            self._calc_magnifications(bad=True)

        model_flux = np.ones(self._dataset.n_epochs) * self.blend_flux
        if self._model.n_sources == 1:
            model_flux += self.source_flux * self._data_magnification
        else:
            for i in range(self._model.n_sources):
                model_flux += self.source_fluxes[i] \
                        * self._data_magnification[i]

        return model_flux

    def get_model_magnitudes(self, **kwargs):
        """
        Arguments:
            **kwargs: see :py:func:`~get_model_fluxes'

        Returns :
            model_mag: *np.ndarray*
                The model magnitude evaluated for each datapoint.
        """
        model_flux = self.get_model_fluxes(**kwargs)
        model_mag = Utils.get_mag_from_flux(model_flux)

        return model_mag

    def scale_fluxes(self, source_flux, blend_flux):
        """
        Rescale the data fluxes to an arbitrary flux scale:
            flux = source_flux_0 * (data.flux - blend_flux) / source_flux
            flux += blend_flux_0
            err_flux = source_flux_0 * data.err_flux / source_flux

        Arguments :
            source_flux: *float*
                Flux of the source in the desired system

            blend_flux: *float*
                Flux of the blend in the desired system

        Returns :
            (flux, err): the fluxes and flux errors of the dataset rescaled
            to the desired system
        """
        if self.model.n_sources > 1:
            raise NotImplementedError(
                'Scaling data to model not implemented for multiple sources.')

        flux = source_flux * (self._dataset.flux - self.blend_flux)
        flux /= self.source_flux
        flux += blend_flux
        err_flux = source_flux * self._dataset.err_flux / self.source_flux

        return (flux, err_flux)

    def get_residuals(
            self, phot_fmt=None, source_flux=None, blend_flux=None, bad=False,
            type=None):
        """
        Calculate the residuals for each datapoint relative to the model.

        Keywords :
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

            type: deprecated, see "phot_fmt".

        Returns :
            residuals: *np.ndarray*
                the residuals for the corresponding dataset.

            errorbars: *np.ndarray*
                the scaled errorbars for each point. For plotting
                errorbars for the residuals.
        """
        if type is not None:
            if type == 'mag':
                warnings.warn(
                    '"mag" returns residuals in the original data flux' +
                    'system. To scale the residuals, use "scaled".')

            warnings.warn(
                'type keyword will be deprecated. Use "phot_fmt" instead.',
                FutureWarning)
            phot_fmt = type

        if bad:
            self._calc_magnifications(bad=True)

        if phot_fmt == 'scaled':
            if self.model.n_sources > 1:
                raise NotImplementedError(
                    'Scaling data to model not implemented for multiple ' +
                    'sources.')

            if source_flux is None or blend_flux is None:
                raise ValueError(
                    'If phot_fmt=scaled, source_flux and blend_flux must' +
                    'also be specified.')

            magnification = self._data_magnification
            model_mag = Utils.get_mag_from_flux(
                blend_flux + source_flux * magnification)
            (flux, err_flux) = self.scale_fluxes(source_flux, blend_flux)
            (mag, err) = Utils.get_mag_and_err_from_flux(flux, err_flux)
            residuals = mag - model_mag
            errorbars = err
        elif phot_fmt == 'mag':
            model_mag = self.get_model_magnitudes()
            residuals = self._dataset.mag - model_mag
            errorbars = self._dataset.err_mag
        elif phot_fmt == 'flux':
            model_flux = self.get_model_fluxes()
            residuals = self._dataset.flux - model_flux
            errorbars = self._dataset.err_flux
        else:
            raise ValueError(
                'phot_fmt must be one of "mag", "flux", or "scaled". Your' +
                'value: {0}'.format(phot_fmt))

        return (residuals, errorbars)

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
                'chi2_gradient() only implemented for single lens models')

        # Implemented for finite source effects?
        if 'rho' in parameters or 't_star' in parameters:
            as_dict = self.model.parameters.as_dict()
            if 'rho' in as_dict or 't_star' in as_dict:
                raise NotImplementedError(
                    'Event.chi2_gradient() is not working ' +
                    'for finite source models yet')

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
        self.calc_chi2_gradient(parameters)
        return self.chi2_gradient


    def calc_chi2_gradient(self, parameters):
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

        # Setup
        gradient = {param: 0 for param in parameters}
        as_dict = self.model.parameters.as_dict()

        # Calculate factor
        # JCY - Everything below here should be refactored into smaller bits.
        factor = self.dataset.flux - self.get_model_fluxes()
        factor *= -2. / self.dataset.err_flux ** 2
        factor *= self.source_flux

        # Get source location
        trajectory = self.model.get_trajectory(self.dataset.time)
        u_2 = trajectory.x ** 2 + trajectory.y ** 2
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
                gradient['t_E'] += np.sum(factor_d_x_d_u * -dt / t_E ** 2)
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
                    as_dict['u_0'] / t_eff ** 2)
        elif 'u_0' not in as_dict:
            t_E = as_dict['t_E'].to(u.day).value
            t_eff = as_dict['t_eff'].to(u.day).value
            if 't_0' in parameters:
                gradient['t_0'] += -sum_d_x_d_u / t_E
            if 't_E' in parameters:
                gradient['t_E'] += (
                                           np.sum(factor_d_x_d_u * dt) -
                                           sum_d_y_d_u * t_eff) / t_E ** 2
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
            det = as_dict['pi_E_N'] ** 2 + as_dict['pi_E_E'] ** 2

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

        self._chi2_gradient = out

        return self._chi2_gradient

    @property
    def chi2_gradient(self):
        """
        Return previously calculated chi^2 gradient (also called Jacobian),
        i.e., :math:`d chi^2/d parameter`. See :py:func:`~get_chi2_gradient()`
        and :py:func:`~calc_chi2_gradient()`.

        Returns :
            gradient: *float* or *np.ndarray*
                chi^2 gradient. Will return None if the chi2 gradient was not
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
        Returns :
            chi2: *float*
                the total chi2 for the fitted dataset. Good points only. See
                :py:obj:`~MulensModel.mulensdata.MulensData.good`.

        If None, you need to run :py:func:`~update()` to execute the
        linear fit and calculate the chi2.
        """
        if self.chi2_per_point is None:
            return None
        else:
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
    def source_flux_ratio(self):
        """
        source_flux_ratio = source_flux_1 / source_flux_0

        Returns :
            source_flux_ratio: *float*
                the ratio of the fitted source fluxes or the value set by
                fix_source_flux_ratio (see :ref:`keywords`).

        If None, you need to run :py:func:`~fit_fluxes()` or
        :py:func:`~update()` to execute the linear fit.
        """
        if self._model.n_sources != 2:
            msg = ("source_flux is defined only for models" +
                   " with TWO sources, you have" +
                   " {0}".format(self._model.n_sources) +
                   " sources.")
            raise NameError(msg)

        if self.fix_source_flux_ratio:
            return self.fix_source_flux_ratio
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
