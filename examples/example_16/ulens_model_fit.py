"""
Class and script for fitting microlensing model using MulensModel.
All the settings are read from a YAML file.
"""
import sys
from os import path
import warnings
import math
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import gridspec, rc, rcParams, rcParamsDefault
from matplotlib.backends.backend_pdf import PdfPages

import_failed = set()
try:
    import yaml
except Exception:
    import_failed.add("yaml")
try:
    import emcee
except Exception:
    import_failed.add("emcee")
try:
    import corner
except Exception:
    import_failed.add("corner")

try:
    import MulensModel as mm
except Exception:
    raise ImportError('\nYou have to install MulensModel first!\n')

__version__ = '0.23.8'


class UlensModelFit(object):
    """
    Class for fitting microlensing model using *MulensModel* package.

    Parameters :
        photometry_files: *list* or *str*
            List of datasets. It can be either a *str* (then just gives
            a name of one file that is read) or a *list*. For *list* each
            element is either a *str* (then just gives the name of the file)
            or a *dict*, which allows more options to be passed, e.g.,

            .. code-block:: python

              [{'file_name': 'data_1.dat', 'phot_fmt': 'mag'}]

            or

            .. code-block:: python

              [{'file_name': 'data_1.dat'}, 'data_2.dat']

            Currently, keyword ``'add_2450000'`` is turned on by default.

        starting_parameters: *dict*
            Starting values of the parameters. Keys of this *dict* are
            microlensing parameters recognized by *MulensModel*. Values
            depend on method used for fitting.

            For EMCEE the values are *str*. First word indicates
            the distribution (allowed: ``gauss``, ``uniform``, and
            ``log-uniform``) and is followed by its parameters.
            For ``uniform`` and ``log-uniform`` these parameters are lower
            and upper limit.
            For ``gauss`` these parameters are mean and sigma of
            the distribution. For example:

            .. code-block:: python

              {
                  't_0': 'gauss 2455703. 0.01',
                  'u_0': 'uniform 0.001 1.',
                  't_E': 'gauss 20. 5.'
              }

        model: *dict*
            Additional settings for *MulensModel.Model*. Accepted keys:

            `'coords'` - event coordinates,

            `'methods'` - methods used for magnification calculation,

            `'methods source 1'` - methods used for magnification calculation
            for the first source in binary source models,

            `'methods source 2'` - methods used for magnification calculation
            for the second source in binary source models,

            `'default method'` - default magnification calculation method,

            `'limb darkening u'` - specifies a *dict* that gives limb
            darkening coefficients in "u" convention, e.g.,
            {'I': 0.4, 'V': 0.5}; note that for plotting the best model we use
            the LD coefficient same as for the first dataset,

            `'parameters'` and `'values'` - used to plot specific model.

        fixed_parameters: *dict*
            Provide parameters that will be kept fixed during the fitting
            process. This option is often used to set parameter reference
            epoch, e.g., ``{'t_0_par': 2456789.}``.

        min_values: *dict*
            Minimum values of parameters that define the prior, e.g.,
            ``{'t_E': 0.}``. Note that the these are only limits of a prior.
            Functional form of priors can be defines in ``fit_constraints``.

        max_values: *dict*
            Maximum values of parameters that define the prior, e.g.,
            ``{'u_0': 1.}``.

        fitting_parameters: *dict*
            Parameters of the fit function. They depend on the method used.

            For EMCEE, the required parameter is ``n_steps``.
            You can also specify ``n_burn`` and ``n_walkers``. The ``n_burn``
            controls the length of burn-in. If not provided, it is assumed to
            be ``0.25*n_steps``. The ``n_walkers`` gives number of parallel
            walkers to be run. If not provided, it is assumed four times
            the number of parameters to be fitted.
            Other options are described below.

            It is possible to export posterior to a .npy file. Just provide
            the file name as ``posterior file`` parameter. You can read this
            file using ``numpy.load()``. You will get an array with a shape of
            (``n_walkers``, ``n_steps-n_burn``, ``n_parameters``). You can
            additionally add option ``posterior file fluxes`` for which
            allowed values are ``all`` and *None* (``null`` in yaml file).
            The value ``all`` means that additionally all source and blending
            fluxes will be saved (``n_parameters`` increases by two times the
            number of datasets).

        fit_constraints: *dict*
            Constraints on model other than minimal and maximal values.

            Currently accepted keys:

            ``'no_negative_blending_flux'`` - reject models with negative
            blending flux if *True*

            ``'negative_blending_flux_sigma_mag'`` - impose a prior that
            disfavours models with negative blending flux using gaussian prior
            for negative values; the value provided should be on the order of
            *20.*

            ``'prior'`` - specifies the priors for quantities. It's also
            a *dict*. Possible key-value pairs:

                ``'t_E': 'Mroz et al. 2017'`` - efficiency-corrected t_E
                distribution from that paper with two modifications: 1) it is
                constant for t_E < 1d, 2) it follows Mao & Paczynski (1996)
                analytical approximation (i.e., slope of -3) for t_E longer
                than probed by Mroz et al. (2017; i.e., 316 d). Note that
                Mroz et al. (2020) studied Galactic bulge.

                ``'t_E': 'Mroz et al. 2020'`` - similar to above but for
                Mroz et al. (2020), where Galactic disc outside bulge region
                was studied. Approximate slopes of 3 and -3 from
                Mao & Paczynski (1996) are used for t_E shorter and longer,
                respectively, than probed by Mroz et al. (2020).

                ``'pi_E_N': gauss mean sigma`` (same for ``'pi_E_E'``) -
                specify gaussian prior for parallax components. Parameters
                *mean* and *sigma* are floats.

            References:
              Mao & Paczynski 1996 -
              https://ui.adsabs.harvard.edu/abs/1996ApJ...473...57M/abstract
              Mroz et al. 2017 -
              https://ui.adsabs.harvard.edu/abs/2017Natur.548..183M/abstract
              Mroz et al. 2020 -
              https://ui.adsabs.harvard.edu/abs/2020ApJS..249...16M/abstract

        plots: *dict*
            Parameters of the plots to be made after the fit. Currently
            allowed keys are ``'triangle'`` and ``'best model'``.
            The values are also dicts and currently accepted keys are
            ``'file'`` (both plots) and ``'time range'`` and for best model
            plot also: ``'magnitude range'``, ``'legend'``, ``'rcParams'``,
            e.g.,

            .. code-block:: python

              {
                  'triangle': {'file': 'my_fit_triangle.png'},
                  'best model':
                      'file': 'my_fit_best.png'
                      'time range': 2456000. 2456300.
                      'magnitude range': 15.123 13.012
                      'legend':
                          'ncol': 2
                          'loc': 'lower center'
                      'rcParams':
                          'font.size': 15
              }

            Note that 'rcParams' allows setting many matplotlib parameters.

        other_output: *dict*
            Parameters for other output. Currently, the only allowed value is
            ``'models': {'file name': NAME_OF_FILE}`` where NAME_OF_FILE is
            a *str* that gives a path to text file to which we will print all
            models and their chi^2. If ``NAME_OF_FILE`` is ``"-"``, then
            the models will be printed to standard output.
    """
    def __init__(
            self, photometry_files, starting_parameters=None, model=None,
            fixed_parameters=None,
            min_values=None, max_values=None, fitting_parameters=None,
            fit_constraints=None, plots=None, other_output=None
            ):
        self._check_MM_version()
        self._photometry_files = photometry_files
        self._starting_parameters = starting_parameters
        self._model_parameters = model
        self._fixed_parameters = fixed_parameters
        self._min_values = min_values
        self._max_values = max_values
        self._fitting_parameters = fitting_parameters
        self._fit_constraints = fit_constraints
        self._plots = plots
        self._other_output = other_output

        self._which_task()
        self._set_default_parameters()
        self._check_imports()

    def _check_MM_version(self):
        """
        Check if MulensModel is new enough
        """
        if int(mm.__version__.split('.')[0]) < 2:
            raise RuntimeError(
                "ulens_model_fit.py requires MulensModel in version "
                "at least 2.0, but you are using " + mm.__version__)

    def _which_task(self):
        """
        Check if input parameters indicate run_fit() or plot_best_model() will
        be run.
        """
        if self._starting_parameters is not None:
            fit = True
        else:
            fit = False
        plot = False

        if self._model_parameters is not None:
            keys = set(self._model_parameters.keys())
            check = keys.intersection({'parameters', 'values'})
            if len(check) == 1:
                raise ValueError(
                    'You have to specify either both or none of ' +
                    'model["parameters"] and model["values"].')
            if len(check) == 2:
                plot = True

        if plot and fit:
            raise ValueError(
                'Too many parameters specified!\nThe starting_parameters ' +
                'indicate you want to fit, but model["parameters"] and ' +
                'model["values"] indicate you want to plot. Please decide')

        if not plot and not fit:
            if self._fixed_parameters is None:
                raise ValueError(
                    'Missing input information. Please specify parameters ' +
                    'to be plotted (model["parameters"] and ' +
                    'model["values"]) or starting_parameters to be fit.')
            else:
                plot = True

        if fit:
            self._task = 'fit'
        elif plot:
            self._task = 'plot'
            self._check_unnecessary_settings()
        else:
            raise ValueError('internal error')

    def _check_unnecessary_settings(self):
        """
        Make sure that there arent' too many parameters specified
        """
        keys = ['_starting_parameters', '_min_values', '_max_values',
                '_fitting_parameters']
        for key in keys:
            if getattr(self, key) is not None:
                raise ValueError(
                    'In plotting mode you should not provide in __init__: ' +
                    key[1:])

        if self._plots is not None:
            if "triangle" in self._plots:
                raise ValueError(
                    'You cannot provide plots["triangle"] if you ' +
                    "don't fit")

    def _set_default_parameters(self):
        """
        set some default parameters
        """
        if self._task == 'fit':
            self._fit_method = 'emcee'
            self._flat_priors = True  # Are priors only 0 or 1?
            self._return_fluxes = True
            self._best_model_ln_prob = -np.inf
        elif self._task == 'plot':
            pass
        else:
            raise ValueError('internal error - task ' + str(self._task))
        self._print_model = False

    def _check_imports(self):
        """
        check if all the required packages are imported
        """
        required_packages = set()

        if self._task == 'fit':
            if self._fit_method == 'emcee':
                required_packages.add('emcee')
            if self._plots is not None and 'triangle' in self._plots:
                required_packages.add('corner')

        failed = import_failed.intersection(required_packages)

        if len(failed) > 0:
            message = (
                'Some of the required packages could not be imported:\n' +
                " ".join(failed))
            if "corner" in failed:
                message += (
                    "\nFor corner package it's enough that you run:\nwget " +
                    "https://raw.githubusercontent.com/dfm/corner.py/" +
                    "v2.0.0/corner/corner.py")
            raise ImportError(message)

    def run_fit(self):
        """
        Run the fit, print the output, and make the plots.

        This function does not accept any parameters. All the settings
        are passed via __init__().
        """
        if self._task != "fit":
            raise ValueError('wrong settings to run .run_fit()')

        self._check_plots_parameters()
        self._check_model_parameters()
        self._parse_other_output_parameters()
        self._get_datasets()
        self._get_parameters_ordered()
        self._get_parameters_latex()
        self._parse_fitting_parameters()
        self._get_n_walkers()
        self._set_min_max_values()
        self._parse_fit_constraints()
        self._parse_starting_parameters()
        self._check_fixed_parameters()
        self._make_model_and_event()
        self._generate_random_parameters()
        self._setup_fit()
        self._run_fit()
        self._finish_fit()
        self._parse_results()
        self._make_plots()

    def plot_best_model(self):
        """
        Plot the best model.
        """
        # XXX - NOTE how the model is defined

        if self._task != "plot":
            raise ValueError('wrong settings to run .plot_best_model()')

        self._check_plots_parameters()
        self._check_model_parameters()
        self._get_datasets()
        self._check_fixed_parameters()
        self._make_model_and_event()
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
                'Unknown plot types: {:}\n'.format(unknown) +
                'Accepted plot types are: ' + ", ".join(allowed_keys))

        for (key, value) in self._plots.items():
            if value is None:
                self._plots[key] = dict()

        if 'best model' in self._plots:
            self._check_plots_parameters_best_model()

        if 'triangle' in self._plots:
            self._check_plots_parameters_triangle()

        if 'best model' in self._plots and 'triangle' in self._plots:
            file_1 = self._plots['best model'].get('file')
            file_2 = self._plots['triangle'].get('file')
            if file_1 == file_2 and file_1 is not None:
                raise ValueError(
                    'Output files for "best model" and "triangle" plots '
                    'cannot be identical')

    def _check_plots_parameters_best_model(self):
        """
        Check if parameters of best model make sense
        """
        allowed = set(['file', 'time range', 'magnitude range', 'legend',
                       'rcParams', 'second Y scale'])
        unknown = set(self._plots['best model'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "best model": {:}'.format(unknown))

        if 'time range' in self._plots['best model']:
            text = self._plots['best model']['time range'].split()
            if len(text) != 2:
                raise ValueError(
                    "'time range' for 'best model' should specify 2 " +
                    "values (begin and end); got: " +
                    str(self._plots['best model']['time range']))
            t_0 = float(text[0])
            t_1 = float(text[1])
            if t_1 < t_0:
                raise ValueError(
                    "Incorrect 'time range' for 'best model':\n" +
                    text[0] + " " + text[1])
            self._plots['best model']['time range'] = [t_0, t_1]

        if 'magnitude range' in self._plots['best model']:
            text = self._plots['best model']['magnitude range'].split()
            if len(text) != 2:
                raise ValueError(
                    "'magnitude range' for 'best model' should specify 2 " +
                    "values (begin and end); got: " +
                    str(self._plots['best model']['magnitude range']))
            mag_0 = float(text[0])
            mag_1 = float(text[1])
            if mag_1 > mag_0:
                raise ValueError(
                    "Incorrect 'magnitude range' for 'best model':\n" +
                    text[0] + " " + text[1])
            self._plots['best model']['magnitude range'] = [mag_0, mag_1]

        for key in ['legend', 'rcParams', 'second Y scale']:
            if key in self._plots['best model']:
                if not isinstance(self._plots['best model'][key], dict):
                    msg = ('The value of {:} (in best model setttings)'
                           'must be a dictionary, but you provided {:}')
                    args = [key, type(self._plots['best model'][key])]
                    raise TypeError(msg.format(*args))

        if 'second Y scale' in self._plots['best model']:
            self._check_plots_parameters_best_model_Y_scale()

    def _check_plots_parameters_best_model_Y_scale(self):
        """
        Check if parameters for second Y scale make sense.
        This function assumes that the second Y scale will be plotted.
        """
        settings = self._plots['best model']['second Y scale']
        allowed = set(['color', 'label', 'labels', 'magnifications'])
        unknown = set(settings.keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "second Y scale" in '
                '"best model": {:}'.format(unknown))
        if not isinstance(settings['magnifications'], list):
            raise TypeError(
                '"best model" -> "second Y scale" -> "magnifications" has to '
                'be a list, not ' + str(type(settings['magnifications'])))
        for value in settings['magnifications']:
            if not isinstance(value, (int, float)):
                raise TypeError(
                    'Wrong value in magnifications: ' + str(value))
        if 'labels' not in settings:
            settings['labels'] = [
                str(x) for x in settings['magnifications']]
        else:
            if not isinstance(settings['labels'], list):
                raise TypeError(
                    '"best model" -> "second Y scale" -> "labels" has to be '
                    'a list, not ' + str(type(settings['labels'])))
            if len(settings['labels']) != len(settings['magnifications']):
                raise ValueError(
                    'In "best model" -> "second Y scale", labels and '
                    'magnifications must be lists of the same length')

    def _check_plots_parameters_triangle(self):
        """
        Check if parameters of triangle plot make sense
        """
        allowed = set(['file'])
        unknown = set(self._plots['triangle'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "triangle": {:}'.format(unknown))

    def _check_model_parameters(self):
        """
        Check parameters of the MulensModel.Model provided by the user
        directly.
        """
        if self._model_parameters is None:
            self._model_parameters = dict()

        allowed = {'coords', 'default method', 'methods',
                   'methods source 1', 'methods source 2',
                   'parameters', 'values', 'limb darkening u'}
        keys = set(self._model_parameters.keys())
        not_allowed = keys - allowed
        if len(not_allowed) > 0:
            raise ValueError(
                'model keyword is a dict with keys not allowed: ' +
                str(not_allowed))
        for key in {'methods', 'methods source 1', 'methods source 2'}:
            if key in self._model_parameters:
                self._model_parameters[key] = self._parse_methods(
                    self._model_parameters[key])
        check = keys.intersection({'parameters', 'values'})
        if len(check) == 1:
            raise ValueError("If you specify 'parameters' and 'values' for " +
                             "'model', then both have to be defined")
        if len(check) == 2:
            self._model_parameters['parameters'] = (
                self._model_parameters['parameters'].split())
            self._model_parameters['values'] = [
                float(x) for x in self._model_parameters['values'].split()]

        if self._starting_parameters is None:
            return  # Below we do checks valid only for fitting.

        condition_1 = ('pi_E_E' in self._starting_parameters)
        condition_2 = ('pi_E_N' in self._starting_parameters)
        if condition_1 or condition_2:
            if 'coords' not in self._model_parameters:
                raise ValueError("Parallax model requires model['coords'].")

    def _parse_methods(self, methods):
        """
        check if odd elements are floats and parse them
        """
        if isinstance(methods, str):
            _enumerate = enumerate(methods.split())
        elif isinstance(methods, list):
            _enumerate = enumerate(methods)
        else:
            raise TypeError(
                'Wrong type of settings specifying methods used to calculate '
                'magnification ("list" or "str" expected): ' +
                str(type(methods)))

        try:
            out = [float(x) if i % 2 == 0 else x for (i, x) in _enumerate]
        except ValueError:
            raise ValueError(
                "Error in parsing floats in methods:\n" + methods)

        if len(out) < 3 or len(out) % 2 != 1:
            raise ValueError(
                "Error in parsing methods:\n" + methods)

        return out

    def _parse_other_output_parameters(self):
        """
        parse information on other output
        """
        if self._other_output is None:
            return

        for (key, value) in self._other_output.items():
            if key == 'models':
                if not isinstance(value, dict):
                    raise ValueError('models value should also be *dict*, ' +
                                     'got ' + str(type(value)))
                for (key2, value2) in value.items():
                    if key2 == 'file name':
                        self._print_model = True
                        self._print_model_i = 0
                        if value2 == '-':
                            self._print_model_file = sys.stdout
                        else:
                            try:
                                self._print_model_file = open(value2, 'w')
                            except Exception:
                                raise ValueError(
                                    'Error while opening file ' + str(value2))
                    else:
                        raise KeyError("Unrecognized key: " + str(key) +
                                       "\nExpected keys: 'file name'.")
            else:
                raise ValueError('Unrecognized key: ' + str(key) + "\n" +
                                 "Expected keys: models")

    def _get_datasets(self):
        """
        construct a list of MulensModel.MulensData objects
        """
        kwargs = {'add_2450000': True}
        if isinstance(self._photometry_files, str):
            self._photometry_files = [self._photometry_files]
        elif not isinstance(self._photometry_files, list):
            raise TypeError(
                'photometry_files should be a list or a str, but you '
                'provided ' + str(type(self._photometry_files)))
        files = [f if isinstance(f, dict) else {'file_name': f}
                 for f in self._photometry_files]
        self._datasets = []
        for file_ in files:
            try:
                dataset = mm.MulensData(**{**kwargs, **file_})
            except FileNotFoundError:
                raise FileNotFoundError(
                    'Provided file path does not exist: ' +
                    str(file_['file_name']))
            except Exception:
                print('Something went wrong while reading file ' +
                      str(file_['file_name']), file=sys.stderr)
                raise
            self._datasets.append(dataset)

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
            t_0='\\Delta t_0', u_0='u_0',
            t_0_1='\\Delta t_{0,1}', u_0_1='u_{0,1}',
            t_0_2='\\Delta t_{0,2}', u_0_2='u_{0,2}', t_E='t_{\\rm E}',
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

        required = ['n_steps']
        strings = ['posterior file', 'posterior file fluxes']
        allowed = ['n_walkers', 'n_burn'] + strings

        full = required + allowed

        for required_ in required:
            if required_ not in settings:
                raise ValueError('EMCEE method requires fitting parameter: ' +
                                 required_)

        if len(set(settings.keys()) - set(full)) > 0:
            raise ValueError('Unexpected fitting parameters: ' +
                             str(set(settings.keys()) - set(full)))

        for (p, value) in settings.items():
            if not isinstance(value, int) and p not in strings:
                raise ValueError(
                    'Fitting parameter ' + p + ' requires int value; got: ' +
                    str(value) + ' ' + str(type(value)))

        if 'n_burn' in settings:
            if settings['n_burn'] >= settings['n_steps']:
                raise ValueError('You cannot set n_burn >= n_steps.')
        else:
            settings['n_burn'] = int(0.25*self._fitting_parameters['n_steps'])

        if 'posterior file' not in settings:
            self._posterior_file_name = None
            if 'posterior file fluxes' in settings:
                raise ValueError('You cannot set "posterior file fluxes" ' +
                                 'without setting "posterior file"')
        else:
            name = settings['posterior file']
            if not isinstance(name, str):
                raise ValueError('"posterior file" must be string, got: ' +
                                 str(type(name)))
            if name[-4:] != '.npy':
                raise ValueError('"posterior file" must end in ".npy", ' +
                                 'got: ' + name)
            if path.exists(name):
                if path.isfile(name):
                    msg = "Exisiting file " + name + " will be overwritten"
                    warnings.warn(msg)
                else:
                    raise ValueError("The path provided for posterior (" +
                                     name + ") exsists and is a directory")
            self._posterior_file_name = name[:-4]
            self._posterior_file_fluxes = None

        if 'posterior file fluxes' in settings:
            fluxes_allowed = ['all', None]
            if settings['posterior file fluxes'] not in fluxes_allowed:
                raise ValueError('Unrecognized "posterior file fluxes": ' +
                                 settings['posterior file fluxes'])
            self._posterior_file_fluxes = settings['posterior file fluxes']

    def _get_n_walkers(self):
        """
        guess how many walkers and starting values there are
        """
        self._n_walkers = None

        if self._fit_method == 'emcee':
            if 'n_walkers' in self._fitting_parameters:
                self._n_walkers = self._fitting_parameters['n_walkers']
            else:
                self._n_walkers = 4 * len(self._starting_parameters)
        else:
            raise ValueError('internal bug')

        if self._n_walkers is None:
            raise ValueError(
                "Couldn't guess the number of walkers based on " +
                "method: {:}\n".format(self._fit_method) +
                "fitting_parameters: " + str(self._fitting_parameters))

    def _set_min_max_values(self):
        """
        Parse min and max values of parameters so that they're properly
        indexed.
        """
        if self._min_values is None:
            self._min_values = []
        if self._max_values is None:
            self._max_values = []

        for key in self._min_values:
            if key in self._max_values:
                if self._min_values[key] >= self._max_values[key]:
                    fmt = (
                        "This doesn't make sense: for {:} the lower limit " +
                        "is larger than the upper limit: {:} vs {:}")
                    raise ValueError(fmt.format(
                        key, self._min_values[key], self._max_values[key]))

        self._min_values_indexed = self._parse_min_max_values_single(
            self._min_values)
        self._max_values_indexed = self._parse_min_max_values_single(
            self._max_values)

    def _parse_min_max_values_single(self, limits):
        """
        change dict that has str as key to index as key
        """
        out = dict()
        if len(limits) == 0:
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
        self._prior_t_E = None
        self._priors = None

        if self._fit_constraints is None:
            self._fit_constraints = {"no_negative_blending_flux": False}
            return

        if isinstance(self._fit_constraints, list):
            raise TypeError(
                "In version 0.5.0 we've changed type of 'fit_constraints' " +
                "from list to dict. Please correct you input and re-run " +
                "the code. Most probably what you need is:\n" +
                "fit_constraints = {'no_negative_blending_flux': True}")

        allowed_keys_flux = {
            "no_negative_blending_flux", "negative_blending_flux_sigma_mag"}
        allowed_keys = {*allowed_keys_flux, "prior"}
        used_keys = set(self._fit_constraints.keys())
        if len(used_keys - allowed_keys) > 0:
            raise ValueError(
                'unrecognized constraint: {:}'.format(forbidden))
        if len(used_keys.intersection(allowed_keys_flux)) == 2:
            raise ValueError(
                'you cannot specify both no_negative_blending_flux and ' +
                'negative_blending_flux_sigma_mag')
        if "no_negative_blending_flux" not in self._fit_constraints:
            self._fit_constraints["no_negative_blending_flux"] = False

        key = "negative_blending_flux_sigma_mag"
        if key in used_keys:
            self._fit_constraints[key] = mm.Utils.get_flux_from_mag(
                self._fit_constraints[key])

        if 'prior' in self._fit_constraints:
            self._parse_fit_constraints_prior()

    def _parse_fit_constraints_prior(self):
        """
        Check if priors in fit constraint are correctly defined.
        """
        priors = dict()
        for (key, value) in self._fit_constraints['prior'].items():
            if key == 't_E':
                if value == "Mroz et al. 2017":
                    self._prior_t_E = 'Mroz+17'
                elif value == "Mroz et al. 2020":
                    self._prior_t_E = 'Mroz+20'
                else:
                    raise ValueError("Unrecognized t_E prior: " + value)
                self._read_prior_t_E_data()
            elif key in ['pi_E_E', 'pi_E_N']:
                words = value.split()
                if len(words) != 3 or words[0] != 'gauss':
                    msg = "Something went wrong in parsing prior for "
                    msg += "{:}: {:}"
                    raise ValueError(msg.format(key, value))
                try:
                    settings = [words[0], float(words[1]), float(words[2])]
                except Exception:
                    raise ValueError('error in parsing: ' + words[1] + " " +
                                     words[2])
                if settings[2] < 0.:
                    raise ValueError('sigma cannot be negative: ' + words[2])
                priors[key] = settings
            else:
                raise KeyError(
                    "Unrecognized key in fit_constraints/prior: " + key)
            self._flat_priors = False

        if len(priors) > 0:
            self._priors = priors

    def _read_prior_t_E_data(self):
        """
        read data that specify t_E prior and parse them appropriately
        """
        self._prior_t_E_data = dict()

        if self._prior_t_E == 'Mroz+17':
            x = np.array([
                -0.93, -0.79, -0.65, -0.51, -0.37, -0.23, -0.09, 0.05, 0.19,
                0.33, 0.47, 0.61, 0.75, 0.89, 1.03, 1.17, 1.31, 1.45, 1.59,
                1.73, 1.87, 2.01, 2.15, 2.29, 2.43])
            y = np.array([
                299.40, 245.60, 358.50, 116.96, 0.00, 47.78, 85.10, 90.50,
                315.37, 501.77, 898.26, 1559.68, 2381.46, 2849.11, 3405.00,
                3431.30, 3611.76, 3038.06, 2170.67, 1680.38, 814.70, 444.06,
                254.89, 114.19, 52.14])
            dx = x[1] - x[0]
            x_min = 0.
            x_max = x[-1] + 0.5 * dx
            mask = (x > x_min-dx)  # We need one more point for extrapolation.
            function = interp1d(x[mask], np.log(y[mask]),
                                kind='cubic', fill_value="extrapolate")
            self._prior_t_E_data['x_min'] = x_min
            self._prior_t_E_data['x_max'] = x_max
            self._prior_t_E_data['y_min'] = function(x_min)
            self._prior_t_E_data['y_max'] = function(x_max)
            self._prior_t_E_data['function'] = function
        elif self._prior_t_E == 'Mroz+20':
# XXX - TO DO:
# - documentation
# - smooth the input data from M+20 and note that
            x = np.array([
                0.74, 0.88, 1.01, 1.15, 1.28, 1.42, 1.55, 1.69, 1.82, 1.96,
                2.09, 2.23, 2.36, 2.50, 2.63])
            y = np.array([
                82.04, 94.98, 167.76, 507.81, 402.08, 681.61, 1157.51,
                1132.80, 668.12, 412.20, 236.14, 335.34, 74.88, 52.64, 97.78])
            dx = (x[1] - x[0]) / 2.
            x_min = x[0] - dx
            x_max = x[-1] + dx
            function = interp1d(x, np.log(y),
                                kind='cubic', fill_value="extrapolate")
            self._prior_t_E_data['x_min'] = x_min
            self._prior_t_E_data['x_max'] = x_max
            self._prior_t_E_data['y_min'] = function(x_min)
            self._prior_t_E_data['y_max'] = function(x_max)
            self._prior_t_E_data['function'] = function
        else:
            raise ValueError('unexpected internal error')

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
                    'starting parameter: ' + words[0] + ' is not recognized.' +
                    'Allowed parameters: ' + str(accepted_types))
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

        if self._task == 'plot':
            return

        repeated = set(self._fit_parameters).intersection(fixed)
        if len(repeated) > 0:
            raise ValueError(
                'Some parameters are both fitted and fixed: ' +
                '{:}'.format(repeated))

    def _make_model_and_event(self):
        """
        Set internal MulensModel instances: Model and Event
        """
        parameters = self._get_example_parameters()

        kwargs = dict()
        if 'coords' in self._model_parameters:
            kwargs['coords'] = self._model_parameters['coords']

        try:
            self._model = mm.Model(parameters, **kwargs)
        except Exception:
            print("Initializer of MulensModel.Model failed.")
            print("Parameters passed: {:}".format(parameters))
            raise
        self._models_satellite = []
        for dataset in self._datasets:
            if dataset.ephemerides_file is None:
                continue
            model = mm.Model(
                parameters, ephemerides_file=dataset.ephemerides_file,
                **kwargs)
            self._models_satellite.append(model)
        key = 'limb darkening u'
        for model in [self._model] + self._models_satellite:
            if key in self._model_parameters:
                for (band, u_value) in self._model_parameters[key].items():
                    model.set_limb_coeff_u(band, u_value)
            if 'default method' in self._model_parameters:
                model.set_default_magnification_method(
                    self._model_parameters['default method'])
            if 'methods' in self._model_parameters:
                model.set_magnification_methods(
                    self._model_parameters['methods'])
            if 'methods source 1' in self._model_parameters:
                self._model.set_magnification_methods(
                    self._model_parameters['methods source 1'], 1)
            if 'methods source 2' in self._model_parameters:
                self._model.set_magnification_methods(
                    self._model_parameters['methods source 2'], 2)

        self._event = mm.Event(self._datasets, self._model)
        self._event.sum_function = 'numpy.sum'

        self._get_n_fluxes()

    def _get_n_fluxes(self):
        """
        find out how many flux parameters there are
        """
        try:
            self._event.get_chi2()
        except ValueError:
            if 'x_caustic_in' in self._model.parameters.parameters:
                self._model.parameters.x_caustic_in = (
                    self._model.parameters.x_caustic_out + 0.01)
                self._event.get_chi2()
            else:
                raise

        n = 0
        for (i, dataset) in enumerate(self._datasets):
            k = len(self._event.fits[i].source_fluxes) + 1
            # Plus 1 is for blending.
            if i == 0:
                self._n_fluxes_per_dataset = k
            elif k != self._n_fluxes_per_dataset:
                raise ValueError(
                    'Strange internal error with number of source fluxes: ' +
                    "{:} {:} {:}".format(i, k, self._n_fluxes_per_dataset))
            n += k

        self._n_fluxes = n

    def _get_example_parameters(self):
        """
        Generate parameters *dict* according to provided starting and fixed
        parameters.
        """
        parameters = dict()
        # XXX it should be either that we have parameters/values or
        # _starting_parameters
        if self._starting_parameters is None:
            keys = self._model_parameters['parameters']
            values = self._model_parameters['values']
            for (key, value) in zip(keys, values):
                parameters[key] = value
            # XXX this is some kind of a hack:
            self._best_model_theta = []
            self._fit_parameters = []
        else:
            for (key, value) in self._starting_parameters.items():
                # We treat Cassan08 case differently so that
                # x_caustic_in is different than x_caustic_out.
                if key == "x_caustic_in":
                    if value[0] == 'gauss':
                        parameters[key] = (
                            value[1] + value[2] * np.random.randn(1)[0])
                    elif value[0] in ['uniform', 'log-uniform']:
                        parameters[key] = 0.25 * value[1] + 0.75 * value[2]
                    else:
                        raise ValueError('internal error: ' + value[0])
                else:
                    if value[0] == 'gauss':
                        parameters[key] = value[1]
                    elif value[0] == 'uniform':
                        parameters[key] = (value[1] + value[2]) / 2.
                    elif value[0] == 'log-uniform':
                        parameters[key] = (value[1] + value[2]) / 2.
                    else:
                        raise ValueError('internal error: ' + value[0])

        if self._fixed_parameters is not None:
            for (key, value) in self._fixed_parameters.items():
                parameters[key] = value

        return parameters

    def _generate_random_parameters(self):
        """
        Generate a number of starting parameters values.
        It is checked if parameters are within the prior.
        """
        max_iteration = 20 * self._n_walkers
        if self._fit_constraints["no_negative_blending_flux"]:
            max_iteration *= 5

        starting = []
        for parameter in self._fit_parameters:
            settings = self._starting_parameters[parameter]
            values = self._get_samples_from_distribution(
                max_iteration, settings)
            starting.append(values)

        starting = np.array(starting).T.tolist()

        self._check_generated_random_parameters(starting)

    def _get_samples_from_distribution(self, n, settings):
        """
        Get n samples from a given distribution (settings[0]).
        The meaning and number of settings[1:] depends on particular
        distribution.
        """
        if settings[0] == 'gauss':
            values = settings[2] * np.random.randn(n) + settings[1]
        elif settings[0] == 'uniform':
            values = np.random.uniform(
                low=settings[1], high=settings[2], size=n)
        elif settings[0] == 'log-uniform':
            beg = math.log(settings[1])
            end = math.log(settings[2])
            values = np.exp(np.random.uniform(beg, end, n))
        else:
            raise ValueError('Unrecognized keyword: ' + settings[0])
        return values

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
                "fit_constraints.\nGot " + str(len(out)))
        self._starting_points = out

    def _ln_prob(self, theta):
        """
        Log probability of the model - combines _ln_prior(), _ln_like(),
        and constraints which include fluxes.

        NOTE: we're using np.log(), i.e., natural logarithms.
        """
        ln_prior = self._ln_prior(theta)
        if not np.isfinite(ln_prior):
            return self._return_ln_prob(-np.inf)

        ln_like = self._ln_like(theta)
        if not np.isfinite(ln_like):
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
                return (value, [0.] * self._n_fluxes)
            else:
                return value
        else:
            if self._return_fluxes:
                if fluxes is None:
                    raise ValueError('Unexpected error!')
                return (value, fluxes)
            else:
                return value

    def _set_model_parameters(self, theta):
        """
        Set microlensing parameters of self._model

        Note that if only plotting functions are called,
        then self._fit_parameters and theta are empty.
        """
        for (parameter, value) in zip(self._fit_parameters, theta):
            setattr(self._model.parameters, parameter, value)

    def _ln_prior(self, theta):
        """
        Check if fitting parameters are within the prior.
        Constraints from self._fit_constraints:
         - on blending flux are NOT applied here,
           but in _run_flux_checks_ln_prior(),
         - on t_E are applied here.

        NOTE: we're using np.log(), i.e., natural logarithms.
        """
        inside = 0.
        outside = -np.inf

        for (index, limit) in self._min_values_indexed.items():
            if theta[index] < limit:
                return outside

        for (index, limit) in self._max_values_indexed.items():
            if theta[index] > limit:
                return outside

        ln_prior = inside

        if self._prior_t_E is not None:
            self._set_model_parameters(theta)
            ln_prior += self._ln_prior_t_E()

        if self._priors is not None:
            self._set_model_parameters(theta)
            for (parameter, prior_settings) in self._priors.items():
                if parameter in ['pi_E_N', 'pi_E_E']:
                    # Other parameters can be added here.
                    value = self._model.parameters.parameters[parameter]
                    ln_prior += self._get_ln_prior_for_1_parameter(
                        value, prior_settings)
                else:
                    raise ValueError('prior not handled: ' + parameter)

        return ln_prior

    def _get_ln_prior_for_1_parameter(self, value, settings):
        """
        Calculate ln(prior) for given value and settings
        """
        if settings[0] == 'gauss':
            sigma = settings[2]
            diff = value - settings[1]
            return -0.5*(diff/sigma)**2 - math.log(math.sqrt(2*np.pi)*sigma)
        else:
            raise ValueError('Case not handelded yet: ' + settings[0])

    def _ln_prior_t_E(self):
        """
        Get log prior for t_E of current model. This function is executed
        if there is t_E prior.
        """
        if self._prior_t_E not in ['Mroz+17', 'Mroz+20']:
            raise ValueError('unexpected internal error ' + self._prior_t_E)

        try:
            x = math.log10(self._model.parameters.t_E)
        except ValueError:
            if 'x_caustic_in' in self._model.parameters.parameters:
                return -np.inf
            else:
                raise

        if x > self._prior_t_E_data['x_max']:
            dy = -3. * math.log(10) * (x - self._prior_t_E_data['x_max'])
            return self._prior_t_E_data['y_max'] + dy
        elif x > self._prior_t_E_data['x_min']:
            return self._prior_t_E_data['function'](x)
        else:
            out = self._prior_t_E_data['y_min'] + 0.
            if self._prior_t_E == 'Mroz+20':
                out += 3. * math.log(10) * (x - self._prior_t_E_data['x_min'])
            return out

    def _ln_like(self, theta):
        """
        likelihood function
        """
        self._set_model_parameters(theta)

        chi2 = self._event.get_chi2()

        if self._print_model:
            self._print_current_model(theta, chi2)

        return -0.5 * chi2

    def _print_current_model(self, theta, chi2):
        """
        print the chi2 and parameters for model provided
        """
        out = "{:.4f}  {:}".format(chi2, " ".join([repr(x) for x in theta]))

        flush = False
        cond_1 = self._print_model_i <= 1000 and self._print_model_i % 10 == 0
        cond_2 = self._print_model_i > 1000 and self._print_model_i % 100 == 0
        if self._print_model_i < 100 or cond_1 or cond_2:
            flush = True

        print(out, file=self._print_model_file, flush=flush)
        self._print_model_i += 1

    def _get_fluxes(self):
        """
        Extract all fluxes and return them in a list.
        """
        fluxes = []
        for (i, dataset) in enumerate(self._datasets):
            fluxes += self._event.fits[i].source_fluxes.tolist()
            fluxes.append(self._event.fits[i].blend_flux)

        return fluxes

    def _run_flux_checks_ln_prior(self, fluxes):
        """
        Run the checks on fluxes - are they in the prior?
        """
        inside = 0.
        outside = -np.inf

        if self._fit_constraints["no_negative_blending_flux"]:
            blend_index = self._n_fluxes_per_dataset - 1
            if fluxes[blend_index] < 0.:
                return outside

        key = "negative_blending_flux_sigma_mag"
        if key in self._fit_constraints:
            blend_index = self._n_fluxes_per_dataset - 1
            if fluxes[blend_index] < 0.:
                sigma = self._fit_constraints[key]
                inside += -0.5 * (fluxes[blend_index] / sigma)**2

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

    def _finish_fit(self):
        """
        Make the things that are necessary after the fit is done.
        Currently it's just closing the file with all models.
        """
        if self._print_model:
            if self._print_model_file is not sys.stdout:
                self._print_model_file.close()
                self._print_model = False

    def _parse_results(self):
        """
        Call the function that prints and saves results
        """
        self._parse_results_EMCEE()
        if self._posterior_file_name is not None:
            self._save_posterior_EMCEE()

    def _parse_results_EMCEE(self):
        """
        Print and save results from EMCEE fitting.

        This version works with EMCEE version 2.X and 3.0.
        """
        n_burn = self._fitting_parameters['n_burn']
        n_fit = len(self._fit_parameters)

        self._samples = self._sampler.chain[:, n_burn:, :].reshape((-1, n_fit))
        print("Fitted parameters:")
        self._parse_results_EMECEE_print(self._samples, self._fit_parameters)
        for name in ['t_0', 't_0_1', 't_0_2']:
            if name in self._fit_parameters:
                index = self._fit_parameters.index(name)
                try:
                    self._samples[:, index] -= int(
                        np.mean(self._samples[:, index]))
                except TypeError:
                    fmt = ("Warning: extremely wide range of posterior t_0: " +
                           "from {:} to {:}")
                    warnings.warn(fmt.format(np.min(self._samples[:, index]),
                                             np.max(self._samples[:, index])))
                    mean = int(np.mean(self._samples[:, index]))
                    self._samples[:, index] = self._samples[:, index] - mean

        if self._return_fluxes:
            try:
                blobs = np.array(self._sampler.blobs)
            except Exception as exception:
                raise ValueError('There was some issue with blobs\n' +
                                 str(exception))
            blob_sampler = np.transpose(blobs, axes=(1, 0, 2))
            blob_samples = blob_sampler[:, n_burn:, :].reshape(
                (-1, self._n_fluxes))
            print("Fitted fluxes (source and blending):")
            if self._n_fluxes_per_dataset == 2:
                s_or_b = ['s', 'b']
            elif self._n_fluxes_per_dataset == 3:
                s_or_b = ['s1', 's2', 'b']
            else:
                raise ValueError(
                    'Internal error: ' + str(self._n_fluxes_per_dataset))
            text = 'flux_{:}_{:}'
            n = self._n_fluxes_per_dataset
            flux_names = [
                text.format(s_or_b[i % n], i // n+1)
                for i in range(self._n_fluxes)
                ]
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
        else:
            self._ln_like(self._best_model_theta)
            print("chi2 : {:.4f}".format(self._event.get_chi2()))
        print(*self._fit_parameters)
        print(*list(self._best_model_theta))
        if self._return_fluxes:
            print("Fluxes:")
            print(*list(self._best_model_fluxes))

    def _save_posterior_EMCEE(self):
        """
        save 3D cube with posterior to a numpy array
        """
        n_burn = self._fitting_parameters.get('n_burn', 0)
        samples = self._sampler.chain[:, n_burn:, :]
        if self._posterior_file_fluxes == 'all':
            blobs = np.array(self._sampler.blobs)
            blobs = np.transpose(blobs, axes=(1, 0, 2))[:, n_burn:, :]
            samples = np.dstack((samples, blobs))
        np.save(self._posterior_file_name, samples)

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
        self._reset_rcParams()

        n_bins = 40

        kwargs = {
            'bins': n_bins, 'labels': self._fit_parameters_latex,
            'show_titles': True, 'quantiles': [0.15866, 0.5, 0.84134],
            'verbose': False, 'top_ticks': False}

        figure = corner.corner(self._samples, **kwargs)

        self._save_figure(self._plots['triangle'].get('file'), figure=figure)

    def _reset_rcParams(self):
        """
        Reset matplotlib rcParams to their defaults
        """
        rcParams.update(rcParamsDefault)

    def _save_figure(self, file_name, figure=None, dpi=None):
        """
        Save figure or display it
        """
        if file_name is None:
            plt.show()
        # elif file_name[-4:].upper() == ".PDF":
        #    pdf = PdfPages(file_name)
        #    if figure is None:
        #        figure = plt.gcf()
        #    pdf.savefig(figure)
        else:
            caller = plt
            if figure is not None:
                caller = figure
            kwargs = dict()
            if dpi is not None:
                kwargs = {'dpi': dpi}
            caller.savefig(file_name, **kwargs)
        plt.close()

    def _best_model_plot(self):
        """
        plot best model and residuals
        """
        dpi = 300

        self._ln_like(self._best_model_theta)  # Sets all parameters to
        # the best model.

        self._reset_rcParams()
        if 'rcParams' in self._plots['best model']:
            for (key, value) in self._plots['best model']['rcParams'].items():
                rcParams[key] = value

        kwargs_all = self._get_kwargs_for_best_model_plot()
        (kwargs_grid, kwargs_model, kwargs, xlim, t_1, t_2) = kwargs_all[:6]
        (kwargs_axes_1, kwargs_axes_2) = kwargs_all[6:]
        (ylim, ylim_residuals) = self._get_ylim_for_best_model_plot(t_1, t_2)

        grid = gridspec.GridSpec(**kwargs_grid)

        axes = plt.subplot(grid[0])
        self._event.plot_data(**kwargs)
        fluxes = self._event.get_ref_fluxes()

        # Plot models below, first ground-based (if needed, hence loop),
        # then satellite ones (if needed).
        for dataset in self._datasets:
            if dataset.ephemerides_file is None:
                self._model.plot_lc(
                    source_flux=fluxes[0], blend_flux=fluxes[1],
                    **kwargs_model)
                break
        for model in self._models_satellite:
            model.parameters.parameters = {**self._model.parameters.parameters}
            model.plot_lc(source_flux=fluxes[0], blend_flux=fluxes[1],
                          **kwargs_model)

        self._plot_legend_for_best_model_plot()
        plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)
        axes.tick_params(**kwargs_axes_1)
        if "second Y scale" in self._plots['best model']:
            self._mark_second_Y_axis_in_best_plot()

        axes = plt.subplot(grid[1])
        self._event.plot_residuals(**kwargs)
        plt.xlim(*xlim)
        plt.ylim(*ylim_residuals)
        axes.tick_params(**kwargs_axes_2)

        self._save_figure(self._plots['best model'].get('file'), dpi=dpi)

    def _get_kwargs_for_best_model_plot(self):
        """
        prepare kwargs/settings for best plot model
        """
        plot_size_ratio = 5
        hspace = 0
        tau = 1.5
        remove_245 = True
        default_model = {'color': 'black', 'linewidth': 1.0, 'zorder': np.inf}

        kwargs_grid = {
            'nrows': 2, 'ncols': 1, 'height_ratios': [plot_size_ratio, 1],
            'hspace': hspace}
        kwargs = {'subtract_2450000': remove_245}

        (t_1, t_2) = self._get_time_limits_for_plot(tau)

        kwargs_model = {
            't_start': t_1, 't_stop': t_2, **default_model, **kwargs}
        if self._model.n_sources != 1:
            kwargs_model['source_flux_ratio'] = self._datasets[0]
        if self._datasets[0].bandpass is not None:
            key = 'limb darkening u'
            if self._datasets[0].bandpass in self._model_parameters[key]:
                u = self._model_parameters[key][self._datasets[0].bandpass]
                kwargs_model['gamma'] = mm.Utils.u_to_gamma(u)

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

    def _get_time_limits_for_plot(self, tau):
        """
        find limits for the best model plot
        """
        if 'time range' in self._plots['best model']:
            t_1 = self._plots['best model']['time range'][0]
            t_2 = self._plots['best model']['time range'][1]
            return (t_1, t_2)

        if self._model.n_sources == 1:
            t_1 = self._model.parameters.t_0 - tau * self._model.parameters.t_E
            t_2 = self._model.parameters.t_0 + tau * self._model.parameters.t_E
        elif self._model.n_sources == 2:
            t_1 = self._model.parameters.t_0_1
            t_2 = self._model.parameters.t_0_2
            if t_1 > t_2:
                (t_1, t_2) = (t_2, t_1)
            t_1 -= tau * self._model.parameters.t_E
            t_2 += tau * self._model.parameters.t_E
        else:
            raise ValueError('internal issue: ' + str(self._model.n_sources))

        return (t_1, t_2)

    def _get_ylim_for_best_model_plot(self, t_1, t_2):
        """
        Find Y axis ranges for plots of data and their residuals.
        Use t_1 and t_2 to limit the data considered.
        """
        padding = 0.05

        y_1 = y_3 = np.inf
        y_2 = y_4 = -np.inf
        (f_source_0, f_blend_0) = self._event.get_ref_fluxes()

        for (i, data) in enumerate(self._datasets):
            mask = (data.time >= t_1) & (data.time <= t_2)
            if np.sum(mask) == 0:
                continue

            (flux, flux_err) = self._event.fits[i].scale_fluxes(
                f_source_0, f_blend_0)
            (y_value, y_err) = mm.Utils.get_mag_and_err_from_flux(
                flux, flux_err)
            mask &= np.logical_not(np.isnan(y_value) | (y_err < 0.))
            y_1 = min(y_1, np.min((y_value - y_err)[mask]))
            y_2 = max(y_2, np.max((y_value + y_err)[mask]))

            (residuals, err_mag) = self._event.fits[i].get_residuals(
                phot_fmt='scaled', source_flux=f_source_0,
                blend_flux=f_blend_0)
            mask_ = np.isfinite(residuals[mask])
            y_3 = min(y_3, np.min((residuals - err_mag)[mask][mask_]))
            y_4 = max(y_4, np.max((residuals + err_mag)[mask][mask_]))

        if y_1 == np.inf:  # There are no data points in the plot.
            return (None, [0.1, -0.1])

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

        if 'magnitude range' in self._plots['best model']:
            ylim = self._plots['best model']['magnitude range']

        return (ylim, ylim_residuals)

    def _plot_legend_for_best_model_plot(self):
        """
        advanced call to plt.legend()
        """
        if len(self._datasets) > 1 or 'legend' in self._plots['best model']:
            if 'legend' not in self._plots['best model']:
                plt.legend()
            else:
                try:
                    plt.legend(**self._plots['best model']['legend'])
                except Exception:
                    print("\npyplot.legend() failed with kwargs:")
                    print(self._plots['best model']['legend'], "\n")
                    raise

    def _mark_second_Y_axis_in_best_plot(self):
        """
        Mark the second (right-hand side) scale for Y axis in
        the best model plot
        """
        settings = self._plots['best model']["second Y scale"]
        magnifications = settings['magnifications']
        color = settings.get("color", "red")
        label = settings.get("label", "magnification")
        labels = settings['labels']

        ylim = plt.ylim()
        flux_min = mm.Utils.get_flux_from_mag(ylim[0])
        flux_max = mm.Utils.get_flux_from_mag(ylim[1])

        (source_flux, blend_flux) = self._event.get_ref_fluxes()
        if self._model.n_sources == 1:
            total_source_flux = source_flux
        else:
            total_source_flux = sum(source_flux)
        flux = total_source_flux * magnifications + blend_flux
        if np.any(flux < 0.):
            mask = (flux > 0.)
            flux = flux[mask]
            labels = [l for (l, m) in zip(labels, mask) if m]
            msg = ("\n\n{:} label/s on the second Y scale will not be shown "
                   "because they correspond to negative flux which cannot "
                   "be translated to magnitudes.")
            warnings.warn(msg.format(np.sum(np.logical_not(mask))))
        A_min = (flux_min - blend_flux) / total_source_flux
        A_max = (flux_max - blend_flux) / total_source_flux

        if (np.min(magnifications) < A_min or np.max(magnifications) > A_max or
                np.any(flux < 0.)):
            msg = ("Provided magnifications for the second (i.e., right-hand "
                   "side) Y-axis scale are from {:} to {:},\nbut the range "
                   "of plotted magnifications is from {:} to {:}, hence, "
                   "the second scale is not plotted")
            args = [min(magnifications), max(magnifications),
                    A_min[0], A_max[0]]
            warnings.warn(msg.format(*args))
            return

        ticks = mm.Utils.get_mag_from_flux(flux)

        ax2 = plt.gca().twinx()
        ax2.set_ylabel(label).set_color(color)
        ax2.spines['right'].set_color(color)
        ax2.set_ylim(ylim[0], ylim[1])
        ax2.tick_params(axis='y', colors=color)
        plt.yticks(ticks, labels, color=color)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')
    if 'yaml' in import_failed:
        raise ImportError('module "yaml" could not be imported :(')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = UlensModelFit(**settings)

    ulens_model_fit.run_fit()
