"""
Class and script for fitting microlensing model using MulensModel.
All the settings are read from a YAML file.
"""
import sys
import numpy as np
import yaml
from os import path
from matplotlib import pyplot as plt
from matplotlib import gridspec

import_failed = set()
try:
    import emcee
except Exception:
    import_failed.add("emcee")
try:
    import corner
except Exception:
    import_failed.add("corner")

import MulensModel as mm


__version__ = '0.3.3'


class UlensModelFit(object):
    """
    Class for fitting microlensing model using *MulensModel* package.

    Parameters :
        photometry_files: *list*
            List of datasets. Currently accepts each entry as *dict*, which
            gives settings passed to *MulensModel.MulensData*, e.g.,

            .. code-block:: python

              [{'file_name': 'data_1.dat', 'phot_fmt': 'mag'}]

            or

            .. code-block:: python

              [{'file_name': 'data_1.dat'}, {'file_name': 'data_2.dat'}]

            Currently, ``'add_2450000'`` is turned on by default.

        starting_parameters: *dict*
            Starting values of the parameters. Keys of this *dict* are
            microlensing parameters recognized by *MulensModel*. Values
            depend on method used for fitting.

            For EMCEE the values are *str*. First word indicates
            the distribution (allowed: ``gauss``, ``uniform``, and
            ``log-uniform``) and is followed by its parameters, e.g.,

            .. code-block:: python

              {
                  't_0': 'gauss 2455703. 0.01',
                  'u_0': 'uniform 0.001 1.',
                  't_E': 'gauss 20. 5.'
              }

        fixed_parameters: *dict*
            Provide parameters that will be kept fixed during the fitting
            process. This option is often used to set parameter reference
            epoch, e.g., ``{'t_0_par': 2456789.}``.

        min_values: *dict*
            Minimum values of parameters that define the prior, e.g.,
            ``{'t_E': 0.}``.

        max_values: *dict*
            Maximum values of parameters that define the prior, e.g.,
            ``{'u_0': 1.}``.

        fitting_parameters: *dict*
            Parameters of the fit function. They depend on the method used.

            For EMCEE, the required parameters are ``n_walkers`` and
            ``n_steps``. Allowed parameter is ``n_burn``.

        fit_constraints: *list*
            List of constraints on model other than minimal and maximal values.

            Currently accepted items:

            ``'no_negative_blending_flux'`` - reject models with negative
            blending flux.

        plots: *dict*
            Parameters of the plots to be made after the fit. Currently
            allowed keys are ``'triangle'`` and ``'best model'``.
            The values are also dicts and currently accept only ``'file'``
            key, e.g.,

            .. code-block:: python

              {
                  'triangle': {'file': 'my_fit_triangle.png'},
                  'best model': {'file': 'my_fit_best.png'}
              }
    """
    def __init__(
            self, photometry_files, starting_parameters=None,
            fixed_parameters=None,
            min_values=None, max_values=None, fitting_parameters=None,
            fit_constraints=None, plots=None,
            ):
        self._photometry_files = photometry_files
        self._starting_parameters = starting_parameters
        self._fixed_parameters = fixed_parameters
        self._min_values = min_values
        self._max_values = max_values
        self._fitting_parameters = fitting_parameters
        self._fit_constraints = fit_constraints
        self._plots = plots

        self._fit_method = 'emcee'
        self._flat_priors = True  # Are priors only 0 or 1?
        self._return_fluxes = True

        self._best_model_ln_prob = -np.inf

        self._check_imports()

    def _check_imports(self):
        """
        check if all the required packages are imported
        """
        required_packages = set()

        if self._fit_method == 'emcee':
            required_packages.add('emcee')
        if self._plots is not None and 'triangle' in self._plots:
            required_packages.add('corner')

        failed = import_failed.intersection(required_packages)

        if len(failed) > 0:
            raise ImportError(
                'Some of the required packages could not be imported:\n' +
                " ".join(failed))

    def run_fit(self):
        """
        Run the fit, print the output, and make the plots.

        This function does not accept any parameters. All the settings
        are passed via __init__().
        """
        self._check_plots_parameters()
        self._get_datasets()
        self._get_parameters_ordered()
        self._get_parameters_latex()
        self._get_n_walkers()
        self._parse_fitting_parameters()
        self._set_min_max_values()
        self._parse_fit_constraints()
        self._parse_starting_parameters()
        self._check_fixed_parameters()
        self._make_model_and_event()
        self._generate_random_parameters()
        self._setup_fit()
        self._run_fit()
        self._parse_results()
        self._make_plots()

    def _check_plots_parameters(self):
        """
        Check if parameters of plots make sense
        """
        allowed_keys = set(['best model', 'triangle'])

        if self._plots is None:
            self._plots = dict()
            return

        unknown = set(self._plots.keys()) - allowed_keys
        if len(unknown) > 0:
            raise ValueError(
                'Unknown plot types: {:}'.format(unknown))

    def _get_datasets(self):
        """
        construct a list of MulensModel.MulensData objects
        """
        kwargs = {'add_2450000': True}
        self._datasets = [
            mm.MulensData(**f, **kwargs) for f in self._photometry_files]

    def _get_parameters_ordered(self):
        """
        Order input parameters in some logical way.
        This is useful to make sure the order of printed parameters
        is always the same.
        """
        all_parameters = (
            't_0 u_0 t_0_1 u_0_1 t_0_2 u_0_2 t_E t_eff rho rho_1 rho_2 ' +
            't_star t_star_1 t_star_2 pi_E_N pi_E_E s q alpha ds_dt ' +
            'dalpha_dt x_caustic_in x_caustic_out t_caustic_in t_caustic_out')
# We do not include t_0_par and t_0_kep because
# these are not fitting parameters.
        all_parameters = all_parameters.split()

        parameters = self._starting_parameters.keys()

        unknown = set(parameters) - set(all_parameters)
        if len(unknown) > 0:
            raise ValueError('Unknown parameters: {:}'.format(unknown))

        indexes = [all_parameters.index(p) for p in parameters]

        self._fit_parameters = [all_parameters[i] for i in indexes]

    def _get_parameters_latex(self):
        """
        change self._fit_parameters into latex parameters
        """
        conversion = dict(
            t_0='\\Delta t_0', u_0='u_0', t_0_1='t_0_1', u_0_1='u_0_1',
            t_0_2='t_0_2', u_0_2='u_0_2', t_E='t_{\\rm E}',
            t_eff='t_{\\rm eff}', rho='\\rho', rho_1='\\rho_1',
            rho_2='\\rho_2', t_star='t_{\\star}', t_star_1='t_{\\star,1}',
            t_star_2='t_{\\star,2}', pi_E_N='\\pi_{{\\rm E},N}',
            pi_E_E='\\pi_{{\\rm E},E}', s='s', q='q', alpha='\\alpha',
            ds_dt='ds/dt', dalpha_dt='d\\alpha/dt',
            x_caustic_in='x_{\\rm caustic,in}',
            x_caustic_out='x_{\\rm caustic,out}',
            t_caustic_in='t_{\\rm caustic,in}',
            t_caustic_out='t_{\\rm caustic,out}')

        self._fit_parameters_latex = [
            ('$' + conversion[key] + '$') for key in self._fit_parameters]

    def _get_n_walkers(self):
        """
        guess how many walkers and starting values there are
        """
        self._n_walkers = None

        if self._fit_method == 'emcee':
            self._n_walkers = self._fitting_parameters.get('n_walkers', None)
        else:
            raise ValueError('internal bug')

        if self._n_walkers is None:
            raise ValueError(
                "Couldn't guess the number of walkers based on " +
                "method: {:}\n".format(self._fit_method) +
                "fitting_parameters: " + str(self._fitting_parameters))

    def _parse_fitting_parameters(self):
        """
        run some checks on self._fitting_parameters to make sure that
        the fit can be run
        """
        if self._fit_method == 'emcee':
            self._parse_fitting_parameters_EMCEE()
        else:
            raise ValueError('internal inconsistency')

    def _parse_fitting_parameters_EMCEE(self):
        """
        make sure EMCEE fitting parameters are properly defined
        """
        settings = self._fitting_parameters

        required = ['n_walkers', 'n_steps']
        allowed = ['n_burn']
        full = required + allowed

        for required_ in required:
            if required_ not in settings:
                raise ValueError('EMCEE method requires fitting parameter: ' +
                                 required_)

        if len(set(settings.keys()) - set(full)) > 0:
            raise ValueError('Unexpected fitting parameters: ' +
                             str(set(settings.keys()) - set(full)))

        for (p, value) in settings.items():
            if not isinstance(value, int):
                raise ValueError(
                    'Fitting parameter ' + p + ' requires int value; got: ' +
                    str(value) + ' ' + str(type(value)))

    def _set_min_max_values(self):
        """
        Parse min and max values of parameters so that they're properly
        indexed.
        """
        self._min_values_indexed = self._parse_min_max_values_single(
            self._min_values)
        self._max_values_indexed = self._parse_min_max_values_single(
            self._max_values)

    def _parse_min_max_values_single(self, limits):
        """
        change dict that has str as key to index as key
        """
        out = dict()
        if limits is None:
            return out

        for (key, value) in limits.items():
            if key not in self._fit_parameters:
                raise ValueError(
                    'Key provided in limits: {:}\n'.format(key) +
                    'is not one of the parameters for fitting: ' +
                    '{:}'.format(self._fit_parameters))
            index = self._fit_parameters.index(key)
            out[index] = value
        return out

    def _parse_fit_constraints(self):
        """
        Parse the fitting constraints that are not simple limits on parameters
        """
        out = {
            "no_negative_blending_flux": False,
            }

        if self._fit_constraints is not None:
            for constraint in self._fit_constraints:
                if constraint == "no_negative_blending_flux":
                    out["no_negative_blending_flux"] = True
                else:
                    raise ValueError(
                        'unrecognized constraint: {:}'.format(constraint))

        self._constraints = out
        self._constraints_any = False
        if len(self._constraints) == 1:
            if not self._constraints["no_negative_blending_flux"]:
                self._constraints_any = True

    def _parse_starting_parameters(self):
        """
        replace self._starting_parameters with dict that has values
        [*str*, *float*, ...]
        and make basic checks
        """
        accepted_types = ['gauss', 'uniform', 'log-uniform']

        out = dict()
        for (key, value) in self._starting_parameters.items():
            words = value.split()
            if words[0] not in accepted_types:
                raise ValueError(
                    'starting value: {:} is not recognized'.format(words[0]))
            if len(words) != 3:
                raise ValueError('Expected 3 parameters, got: ' + str(words))
            floats = []
            for word in words[1:]:
                try:
                    floats.append(float(word))
                except Exception:
                    raise ValueError('Expected float, got {:}'.format(word))
            if words[0] == 'gauss':
                if floats[1] < 0.:
                    raise ValueError(
                        'Sigma cannot be negative, got: ' + str(floats[1]))
            if words[0] in ['uniform', 'log-uniform']:
                if floats[1] < floats[0]:
                    raise ValueError(
                        'For uniform distribution, the second parameters ' +
                        'has to be larger than the first one.\n Got ' +
                        '{:} {:}'.format(floats[0], floats[1]))
            out[key] = [words[0]] + floats

        self._starting_parameters = out

    def _check_fixed_parameters(self):
        """
        Check if fixed_parameters make sense
        """
        if self._fixed_parameters is None:
            return

        all_parameters = (
            't_0 u_0 t_0_1 u_0_1 t_0_2 u_0_2 t_E t_eff rho rho_1 rho_2 ' +
            't_star t_star_1 t_star_2 pi_E_N pi_E_E s q alpha ds_dt ' +
            'dalpha_dt x_caustic_in x_caustic_out ' +
            't_caustic_in t_caustic_out ' +
            't_0_par t_0_kep')
        all_parameters = set(all_parameters.split())

        fixed = set(self._fixed_parameters.keys())

        unknown = fixed - all_parameters
        if len(unknown) > 0:
            raise ValueError('Unknown fixed parameters: {:}'.format(unknown))

        repeated = set(self._fit_parameters).intersection(fixed)
        if len(repeated) > 0:
            raise ValueError(
                'Some parameters are both fitted and fixed: ' +
                '{:}'.format(repeated))

    def _make_model_and_event(self):
        """
        Set internal MulensModel instances: Model and Event
        """
        parameters = dict()
        for (key, value) in self._starting_parameters.items():
            if value[0] == 'gauss':
                parameters[key] = value[1]
            elif value[0] == 'uniform':
                parameters[key] = (value[1] + value[2]) / 2.
            elif value[0] == 'log-uniform':
                parameters[key] = (value[1] + value[2]) / 2.

        if self._fixed_parameters is not None:
            for (key, value) in self._fixed_parameters.items():
                parameters[key] = value

        try:
            self._model = mm.Model(parameters)
        except Exception:
            print("Initializer of MulensModel.Model failed.")
            print("Parameters passed: {:}".format(parameters))
            raise

        self._event = mm.Event(self._datasets, self._model)
        self._event.sum_function = 'numpy.sum'

    def _generate_random_parameters(self):
        """
        Generate a number of starting parameters values.
        It is checked if parameters are within the prior.
        """
        max_iteration = 20 * self._n_walkers

        starting = []
        for parameter in self._fit_parameters:
            settings = self._starting_parameters[parameter]
            if settings[0] == 'gauss':
                values = settings[2] * np.random.randn(max_iteration)
                values += settings[1]
            elif settings[0] == 'uniform':
                values = np.random.uniform(
                    low=settings[1], high=settings[2], size=max_iteration)
            elif settings[0] == 'log-uniform':
                beg = np.log(settings[1])
                end = np.log(settings[2])
                values = np.exp(np.random.uniform(beg, end, max_iteration))
            else:
                raise ValueError('Unrecognized keyword: ' + settings[0])
            starting.append(values)

        starting = np.array(starting).T.tolist()

        self._check_generated_random_parameters(starting)

    def _check_generated_random_parameters(self, starting):
        """
        Check if the set of points provided has at least self._n_walkers
        points inside the prior.
        """
        out = []
        for point in starting:
            ln_prob = self._ln_prob(point)
            if self._return_fluxes:
                ln_prob = ln_prob[0]
            if ln_prob > -np.inf:
                out.append(point)
                if len(out) == self._n_walkers:
                    break

        if len(out) < self._n_walkers:
            raise ValueError(
                "Couldn't generate required starting points in a prior. " +
                "Most probably you have to correct at least one of: " +
                "starting_parameters, min_values, max_values, or " +
                "fit_constraints.")
        self._starting_points = out

    def _ln_prob(self, theta):
        """
        Log probability of the model - combines _ln_prior(), _ln_like(),
        and constraints which include fluxes.
        """
        ln_prior = self._ln_prior(theta)
        if not np.isfinite(ln_prior):
            return self._return_ln_prob(-np.inf)

        ln_like = self._ln_like(theta)
        if not np.isfinite(ln_prior):
            return self._return_ln_prob(-np.inf)

        ln_prob = ln_prior + ln_like

        fluxes = self._get_fluxes()

        ln_prior_flux = self._run_flux_checks_ln_prior(fluxes)
        if not np.isfinite(ln_prior_flux):
            return self._return_ln_prob(-np.inf)

        ln_prob += ln_prior_flux

        self._update_best_model(ln_prob, theta, fluxes)

        return self._return_ln_prob(ln_prob, fluxes)

    def _return_ln_prob(self, value, fluxes=None):
        """
        used to parse output of _ln_prob() in order to make that function
        shorter
        """
        if value == -np.inf:
            if self._return_fluxes:
                if self._model.n_sources > 1:
                    raise ValueError(
                        "I'm not sure - for binary source models are there " +
                        "2 fluxes or 1 for each dataset?")
                n_fluxes = (self._model.n_sources + 1) * len(self._datasets)
                return (value, [0.] * n_fluxes)
            else:
                return value
        else:
            if self._return_fluxes:
                if fluxes is None:
                    raise ValueError('Unexpected error!')
                return (value, fluxes)
            else:
                return value

    def _ln_prior(self, theta):
        """
        Check if fitting parameters are within the prior.
        Constraints from self._constraints and NOT applied here.
        """
        inside = 0.
        outside = -np.inf

        for (index, limit) in self._min_values_indexed.items():
            if theta[index] < limit:
                return outside

        for (index, limit) in self._max_values_indexed.items():
            if theta[index] > limit:
                return outside

        return inside

    def _ln_like(self, theta):
        """
        likelihood function
        """
        for (parameter, value) in zip(self._fit_parameters, theta):
            setattr(self._model.parameters, parameter, value)

        chi2 = self._event.get_chi2()

        return -0.5 * chi2

    def _get_fluxes(self):
        """
        Extract all fluxes and return them in a list.
        """
        if self._model.n_sources > 1:
            raise ValueError(
                "I'm not sure - for binary source models are there " +
                "2 fluxes or 1 for each dataset?")

        fluxes = []
        for dataset in self._datasets:
            fluxes.append(self._event.fit.flux_of_sources(dataset)[0])
            fluxes.append(self._event.fit.blending_flux(dataset))

        return fluxes

    def _run_flux_checks_ln_prior(self, fluxes):
        """
        Run the checks on fluxes - are they in the prior
        """
        inside = 0.
        outside = -np.inf

        if self._constraints["no_negative_blending_flux"]:
            blend_index = self._model.n_sources
            if self._model.n_sources > 1:
                raise ValueError(
                    "I'm not sure - for binary source models are there 2 " +
                    "fluxes or 1 for each dataset?")
            if fluxes[blend_index] < 0.:
                return outside

        return inside

    def _update_best_model(self, ln_prob, theta, fluxes):
        """
        Check if the current model is the best one and save information.
        """
        if ln_prob < self._best_model_ln_prob:
            return

        self._best_model_ln_prob = ln_prob
        self._best_model_theta = theta
        self._best_model_fluxes = fluxes

    def _setup_fit(self):
        """
        Setup what is needed for fitting after MulensModel.Event is set up.
        """
        self._setup_fit_EMCEE()

    def _setup_fit_EMCEE(self):
        """
        Setup fit using EMCEE
        """
        n_fit = len(self._fit_parameters)
        self._sampler = emcee.EnsembleSampler(
            self._n_walkers, n_fit, self._ln_prob)

    def _run_fit(self):
        """
        Call the method that does the fit.
        """
        self._run_fit_EMCEE()

    def _run_fit_EMCEE(self):
        """
        Run EMCEE
        """
        self._sampler.run_mcmc(self._starting_points,
                               self._fitting_parameters['n_steps'])

    def _parse_results(self):
        """
        Call the function that prints and saves results
        """
        self._parse_results_EMCEE()

    def _parse_results_EMCEE(self):
        """
        Print and save results from EMCEE fitting.

        This version works with EMCEE version 2.X and 3.0.
        """
        n_burn = self._fitting_parameters.get('n_burn', 0)
        n_fit = len(self._fit_parameters)

        self._samples = self._sampler.chain[:, n_burn:, :].reshape((-1, n_fit))
        print("Fitted parameters:")
        self._parse_results_EMECEE_print(self._samples, self._fit_parameters)
        if 't_0' in self._fit_parameters:
            index = self._fit_parameters.index('t_0')
            self._samples[:, index] -= int(np.mean(self._samples[:, index]))

        if self._return_fluxes:
            try:
                blobs = np.array(self._sampler.blobs)
            except Exception as exception:
                raise ValueError('There was some issue with blobs\n' +
                                 str(exception))
            blob_sampler = np.transpose(blobs, axes=(1, 0, 2))
            n_fluxes = blob_sampler.shape[-1]
            blob_samples = blob_sampler[:, n_burn:, :].reshape((-1, n_fluxes))
            print("Fitted fluxes (source and blending):")
            s_or_b = ['s', 'b']
            text = 'flux_{:}_{:}'
            flux_names = [
                text.format(s_or_b[i % 2], i // 2+1) for i in range(n_fluxes)]
            self._parse_results_EMECEE_print(blob_samples, flux_names)

        self._print_best_model()

    def _parse_results_EMECEE_print(self, data, parameters):
        """
        print mean values and +- 1 sigma for given parameters
        """
        results = np.percentile(data, [15.866, 50, 84.134], axis=0)
        sigma_plus = results[2, :] - results[1, :]
        sigma_minus = results[1, :] - results[0, :]

        for out in zip(parameters, results[1], sigma_plus, sigma_minus):
            format_ = "{:} : {:.5f} +{:.5f} -{:.5f}"
            if out[0] == 'q':
                format_ = "{:} : {:.7f} +{:.7f} -{:.7f}"
            print(format_.format(*out))

    def _print_best_model(self):
        """
        print best model found
        """
        print("Best model:")
        if self._flat_priors:
            print("chi2 : {:.4f}".format(-2. * self._best_model_ln_prob))
        print(*self._fit_parameters)
        print(*list(self._best_model_theta))
        if self._return_fluxes:
            print("Fluxes:")
            print(*list(self._best_model_fluxes))

    def _make_plots(self):
        """
        make plots after fitting: best model, triangle plot, trace plot
        """
        if 'best model' in self._plots:
            self._best_model_plot()
        if 'triangle' in self._plots:
            self._triangle_plot()

    def _triangle_plot(self):
        """
        Make a triangle plot
        """
        n_bins = 40

        kwargs = {
            'bins': n_bins, 'labels': self._fit_parameters_latex,
            'show_titles': True, 'quantiles': [0.15866, 0.5, 0.84134],
            'verbose': False, 'top_ticks': False}

        figure = corner.corner(self._samples, **kwargs)

        if 'file' in self._plots['triangle']:
            figure.savefig(self._plots['triangle']['file'])
        else:
            plt.show()
        plt.close()

    def _best_model_plot(self):
        """
        plot best model and residuals
        """
        dpi = 300
        kwargs_all = self._get_kwargs_for_best_model_plot()
        (kwargs_grid, kwargs_model, kwargs, xlim, t_1, t_2) = kwargs_all[:6]
        (kwargs_axes_1, kwargs_axes_2) = kwargs_all[6:]

        (ylim, ylim_residuals) = self._get_ylim_for_best_model_plot(t_1, t_2)

        grid = gridspec.GridSpec(**kwargs_grid)

        axes = plt.subplot(grid[0])
        self._event.plot_data(**kwargs)
        self._event.plot_model(**kwargs_model)
        if len(self._datasets) > 1:
            plt.legend()
        plt.xlim(*xlim)
        plt.ylim(*ylim)
        axes.tick_params(**kwargs_axes_1)

        axes = plt.subplot(grid[1])
        self._event.plot_residuals(**kwargs)
        plt.xlim(*xlim)
        plt.ylim(*ylim_residuals)
        axes.tick_params(**kwargs_axes_2)

        if 'file' in self._plots['best model']:
            plt.savefig(self._plots['best model']['file'], dpi=dpi)
        else:
            plt.show()
        plt.close()

    def _get_kwargs_for_best_model_plot(self):
        """
        prepare kwargs/settings for best plot model
        """
        plot_size_ratio = 5
        hspace = 0
        tau = 1.5
        remove_245 = True
        model_color = 'black'

        kwargs_grid = {
            'nrows': 2, 'ncols': 1, 'height_ratios': [plot_size_ratio, 1],
            'hspace': hspace}
        kwargs = {'subtract_2450000': remove_245}

        t_1 = self._model.parameters.t_0 - tau * self._model.parameters.t_E
        t_2 = self._model.parameters.t_0 + tau * self._model.parameters.t_E

        kwargs_model = {
            'color': model_color, 't_start': t_1, 't_stop': t_2, **kwargs}
        if kwargs['subtract_2450000']:
            xlim = [t_1-2450000., t_2-2450000.]
        else:
            xlim = [t_1, t_2]

        kwargs_axes_1 = dict(
            axis='both', direction='in', bottom=True, top=True, left=True,
            right=True, labelbottom=False)
        kwargs_axes_2 = {**kwargs_axes_1, 'labelbottom': True}

        return (kwargs_grid, kwargs_model, kwargs, xlim, t_1, t_2,
                kwargs_axes_1, kwargs_axes_2)

    def _get_ylim_for_best_model_plot(self, t_1, t_2):
        """
        Find Y axis ranges for plots of data and their residuals.
        Use t_1 and t_2 to limit the data considered.
        """
        padding = 0.05
        phot_fmt = 'mag'

        y_1 = y_3 = np.inf
        y_2 = y_4 = -np.inf

        for data in self._datasets:
            mask = (data.time >= t_1) & (data.time <= t_2)
            err_mag = data.err_mag[mask]
            y_1 = min(y_1, np.min(data.mag[mask] - err_mag))
            y_2 = max(y_2, np.max(data.mag[mask] + err_mag))

            residuals = self._model.get_residuals(
                data=data, type=phot_fmt)[0][0][mask]
            y_3 = min(y_3, np.min(residuals - err_mag))
            y_4 = max(y_4, np.max(residuals + err_mag))

        dy = padding * (y_2 - y_1)
        dres = padding * (y_4 - y_3)
        ylim = [y_2 + dy, y_1 - dy]
        ylim_r = [y_4 + dres, y_3 - dres]

        # Block below is the same in MulensModel.Model.plot_residuals() in
        # its version 1.15.6.
        ylim_r_max = np.max(np.abs(ylim_r))
        if ylim_r_max > 1.:
            ylim_r_max = 0.5
        ylim_residuals = [ylim_r_max, -ylim_r_max]

        return (ylim, ylim_residuals)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')

    input_file = sys.argv[1]
    input_file_root = path.splitext(input_file)[0]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = UlensModelFit(**settings)

    ulens_model_fit.run_fit()
