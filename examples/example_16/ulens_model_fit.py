"""
Class and script for fitting microlensing model using MulensModel.
All the settings are read from a YAML file.
"""
import sys
from os import path, sep
import tempfile
import shutil
import warnings
import math
import numpy as np
import shlex
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from matplotlib import gridspec, rcParams, rcParamsDefault
# from matplotlib.backends.backend_pdf import PdfPages

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
    from pymultinest.run import run as mn_run
    from pymultinest.analyse import Analyzer
except Exception:
    import_failed.add("pymultinest")
try:
    import ultranest
except Exception:
    import_failed.add("ultranest")
try:
    import plotly.graph_objects as go
except Exception:
    import_failed.add("plotly")
try:
    import MulensModel as mm
except Exception:
    raise ImportError('\nYou have to install MulensModel first!\n')


__version__ = '0.41.0'


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
            Except standard parameters of MulensData, one can additionally
            pass
            ``'scale_errorbars': {'factor': kappa, 'minimum': epsilon}``
            to scale uncertainties.

        fit_method: *str*
            Method of fitting. Currently accepted values are ``EMCEE``,
            ``MultiNest``, and ``UltraNest``. If not provided, the script
            will guess it based on other parameters: ``EMCEE`` is selected
            if ``starting_parameters`` are provided. If ``prior_limits`` are
            provided, either ``MultiNest`` or ``UltraNest`` will be selected
            depending on the ``fitting_parameters``.

            Webpage of each method:
            - EMCEE: https://emcee.readthedocs.io/en/stable/
            - MultiNest: https://johannesbuchner.github.io/PyMultiNest/
            - UltraNest: https://johannesbuchner.github.io/UltraNest/

        starting_parameters: *dict*
            Starting values of the parameters.
            It also indicates the EMCEE fitting mode.
            There are two possibilities for the information provided.

            First, you can provide a name of file with sets of
            parameters and names of parameters. For example:

            .. code-block:: python

              {
                  'file': 'STARTING_POINTS.txt',
                  'parameters': 't_0 u_0 t_E'
              }

            In that case the text file has three columns and
            at least 2 * 3 = 6 rows with values.

            Second, you can provide distribution that defines distributions
            from which values of parameters will be drawn.
            Keys of this *dict* are microlensing parameters recognized by
            *MulensModel* and values are *str*. First word indicates
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

        prior_limits: *dict*
            Upper and lower limits of parameters.
            It only applies to pyMultiNest and UltraNest fitting.

            Keys are MulensModel parameters and values are lists of two floats
            each (alternatively a string giving 2 floats can be provided - see
            example below). Currently, no informative priors are allowed for
            pyMultiNest and UltraNest fitting. Example input:

            .. code-block:: python

              {
                  't_0': [2455379.4, 2455379.76]
                  'u_0': [0.46, 0.65]
                  't_E': "16. 19.6"
              }

        model: *dict*
            Additional settings for *MulensModel.Model*. Accepted keys:

            ``'coords'`` - event coordinates,

            ``'methods'`` - methods used for magnification calculation,

            ``'methods source 1'`` - methods used for magnification calculation
            for the first source in binary source models,

            ``'methods source 2'`` - methods used for magnification calculation
            for the second source in binary source models,

            ``'default method'`` - default magnification calculation method,

            ``'methods parameters'`` - dict of dicts that add more parameters
            that are passed to methods calculating magnification; typical call:
            ``'VBBL': {'accuracy': 0.01}``

            ``'limb darkening u'`` - specifies a *dict* that gives limb
            darkening coefficients in "u" convention, e.g.,
            {'I': 0.4, 'V': 0.5}; note that for plotting the best model we use
            the LD coefficient same as for the first dataset,

            ``'parameters'`` and ``'values'`` - used to plot specific model.

        fixed_parameters: *dict*
            Provide parameters that will be kept fixed during the fitting
            process. This option is often used to set parameter reference
            epoch, e.g., ``{'t_0_par': 2456789.}``.

        min_values: *dict*
            Minimum values of parameters that define the prior, e.g.,
            ``{'t_E': 0.}``. Note that the these are only limits of a prior.
            Functional form of priors can be defined in ``fit_constraints``.
            It works only for EMCEE fitting.

        max_values: *dict*
            Maximum values of parameters that define the prior, e.g.,
            ``{'u_0': 1.}``.
            It works only for EMCEE fitting.

        fitting_parameters: *dict*
            Parameters of the fit function. They depend on the method used -
            we discuss EMCEE, pyMultiNest and UltraNest below.

            First - EMCEE. The required parameter is ``n_steps``.
            You can also specify ``n_burn`` and ``n_walkers``. The ``n_burn``
            controls the length of burn-in. If not provided, it is assumed to
            be ``0.25*n_steps``. The ``n_walkers`` gives number of parallel
            walkers to be run. If not provided, it is assumed four times
            the number of parameters to be fitted.
            Other options are described below.

            The ``progress`` option (*bool* type value; default is *False*)
            controls if a progress bar is shown.

            It is possible to export posterior to a .npy file. Just provide
            the file name as ``posterior file`` parameter. You can read this
            file using ``numpy.load()``. You will get an array with a shape of
            (``n_walkers``, ``n_steps-n_burn``, ``n_parameters``). You can
            additionally add option ``posterior file fluxes`` for which
            allowed values are ``all`` and *None* (``null`` in yaml file).
            The value ``all`` means that additionally all source and blending
            fluxes will be saved (``n_parameters`` increases by two times the
            number of datasets).

            Second - pyMultiNest. There are no required parameters, but a few
            can be provided. Currently accepted ones are:

            ``basename`` (*str*) - common part of the output files produced by
            MultiNest. If you don't provide it, then the output would be
            saved to temporary files and deleted at the end.

            ``multimodal`` (*bool*) - do you want multiple modes in
            the posterior to be detected and reported separately?

            ``n_live_points`` (*int*) - number of live points, default value
            is 400. Also valid for UltraNest.

            ``sampling efficiency`` (*float*) - requested sampling efficiency.
            MultiNest documentation suggests 0.8 (default value) for parameter
            estimation and 0.3 for evidence evaluation.

            ``evidence tolerance`` (*float*) - requested tolerance of ln(Z)
            calculation; default is 0.5 and should work well in most cases.

            Third - UltraNest. There are no required parameters, but a few
            can be provided. Currently accepted ones are:

            ``log directory`` (*str*) - where to store output files. If given,
            there is a check if directory exists. If not given, no outputs
            are saved.

            ``derived parameter names`` (*str*) - names of additional derived
            parameters created by transform. In microlensing, they are usually
            the source(s) and blending fluxes. If not given, they are ignored
            in the transform function.

            ``show_status`` (*bool*) - whether to show integration progress
            as a status line or not. Default is *True*.

            ``min_num_live_points`` (*int*) - minimum number of live points
            throughout the run. Default value is 400.

            ``dlogz`` (*float*) - Target evidence uncertainty, in order to
            obtain a logz error below a threshold. Default value is 0.5.
            It can be increased to allow `min_num_live_points` values below:
            sqrt(iterations) / dlogz = sqrt(1000) / 0.5 ~ 64.

            ``frac_remain`` (*float*) - Integrate until this fraction of the
            integral is left in the remainder. Numbers smaller than 0.01
            ensure that peaks are discovered, higher numbers can be set if
            the posterior is simple. Default value is 0.01.

            ``max_num_improvement_loops`` (*int*) - Limit the number of times
            the algorithm is repeated to improve. Default value is -1.

        fit_constraints: *dict*
            Constraints on model other than minimal and maximal values.

            Currently accepted keys:

            ``'no_negative_blending_flux'`` - reject models with negative
            blending flux if *True*

            ``'negative_blending_flux_sigma_mag'`` - impose a prior that
            disfavors models with negative blending flux using gaussian prior
            for negative values; the value provided should be on the order of
            *20.*

            ``'color'`` - specify gaussian prior for colors of the sources.
            Parameters:
                *mean* and *sigma* are floats in magnitudes, *dataset_label*
                are str defined in MulensData.plot_properties['label']

            ``'color source 1'`` - specify gaussian prior for color of
            the primary source in binary source model.
            Parameters:
                *mean* and *sigma* are floats in magnitudes, *dataset_label*
                are str defined in MulensData.plot_properties['label']

            ``'color source 2'`` - specify gaussian prior for color of
            the secondary source in binary source model.
            Parameters:
                *mean* and *sigma* are floats in magnitudes, *dataset_label*
                are str defined in MulensData.plot_properties['label']

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

            ``'posterior parsing'`` - additional settings that allow
            modyfing posterior after it's calculated. Possile values:

                ``'abs': [...]`` - calculate absolute values for parameters
                from given list. It's useful for e.g. ``'u_0'`` for
                free-floating planet events.

        plots: *dict*
            Parameters of the plots to be made after the fit. Currently
            allowed keys are ``triangle``, ``trace`` (only EMCEE fitting),
            ``trajectory``, and ``best model``.
            The values are also dicts and currently accepted keys are:
            1) for ``best model``:
            ``'file'``,``'interactive' ``'time range'``, ``'magnitude range'``,
            ``'legend'``,`and ``'rcParams'``,
            2) for ``triangle`` and ``trace``:
            ``'file'`` and ``'shift t_0'`` (*bool*, *True* is default)
            3) for ``trajectory``:
            ``'file'`` and ``'time range'`` (if not provided, then values
            from ``best model`` will be used)
            e.g.:

            .. code-block:: python

              {
                  'triangle':
                      'file': 'my_fit_triangle.png'
                  'trace':
                      'file': 'my_fit_trace_plot.png'
                      'shift t_0': False
                  'trajectory':
                      'file': 'my_trajectory.png'
                      'time range': 2456050. 2456300.
                  'best model':
                      'file': 'my_fit_best.png'
                      'interactive' : 'my_fit_best.html'
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
            Parameters for other output. Allowed value are:

            ``'models': {'file name': NAME_OF_FILE}`` where NAME_OF_FILE is
            a *str* that gives a path to text file to which we will print all
            models and their chi^2. If ``NAME_OF_FILE`` is ``"-"``, then
            the models will be printed to standard output.

            ``'yaml output': {'file name': NAME_OF_FILE}`` where NAME_OF_FILE
            is a *str* that gives a path to YAML-format file to which
            the results will be printed

            ``'residuals': {'files': FILE_1 FILE_2 ...}`` where FILE_X are
            the names of the files to which residuals will be printed.
            These files will have three columns: time, residuals, and
            uncertainties.
    """

    def __init__(
            self, photometry_files,
            starting_parameters=None, prior_limits=None, model=None,
            fixed_parameters=None,
            min_values=None, max_values=None, fitting_parameters=None,
            fit_constraints=None, plots=None, other_output=None,
            fit_method=None
    ):
        self._check_MM_version()
        self._photometry_files = photometry_files
        self._starting_parameters_input = starting_parameters
        self._prior_limits = prior_limits
        self._model_parameters = model
        self._fixed_parameters = fixed_parameters
        self._min_values = min_values
        self._max_values = max_values
        self._fitting_parameters = fitting_parameters
        self._fit_constraints = fit_constraints
        self._plots = plots
        self._other_output = other_output
        self._fit_method = fit_method

        self._which_task()
        self._set_default_parameters()
        if self._task == 'fit':
            if self._fit_method is None:
                self._guess_fitting_method()
            else:
                self._fit_method = self._fit_method.lower()
                self._check_fitting_method()
            self._check_starting_parameters_type()
            self._set_fit_parameters_unsorted()
        self._check_imports()

    def _check_MM_version(self):
        """
        Check if MulensModel is new enough
        """
        # code_version = "{:} and {:}".format(mm.__version__, __version__)
        # print('\nMulensModel and script versions:', code_version, end='\n\n')
        if int(mm.__version__.split('.')[0]) < 2:
            raise RuntimeError(
                "ulens_model_fit.py requires MulensModel in version "
                "at least 2.0, but you are using " + mm.__version__)

    def _which_task(self):
        """
        Check if input parameters indicate run_fit() or plot_best_model() will
        be run.
        """
        if self._starting_parameters_input is not None:
            fit = True
        elif self._prior_limits is not None:
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
            self._check_unnecessary_settings_plot()
        else:
            raise ValueError('internal error')

    def _check_unnecessary_settings_plot(self):
        """
        Make sure that there aren't too many parameters specified for:
        self._task == 'plot'
        """
        keys = ['_starting_parameters_input', '_min_values', '_max_values',
                '_fitting_parameters', '_prior_limits']
        for key in keys:
            if getattr(self, key) is not None:
                raise ValueError(
                    'In plotting mode you should not provide in __init__: ' +
                    key[1:])

        if self._plots is not None:
            if "triangle" in self._plots:
                raise ValueError(
                    'You cannot provide plots["triangle"] if you '
                    "don't fit")
            if "trace" in self._plots:
                raise ValueError(
                    'You cannot provide plots["trace"] if you' " don't fit")

    def _set_default_parameters(self):
        """
        set some default parameters
        """
        if self._task == 'fit':
            self._flat_priors = True  # Are priors only 0 or 1?
            self._return_fluxes = True
            self._best_model_ln_prob = -np.inf
            self._flux_names = None
            self._shift_t_0 = True
        elif self._task == 'plot':
            pass
        else:
            raise ValueError('internal error - task ' + str(self._task))
        self._print_model = False
        self._yaml_results = False
        self._residuals_output = False

        parameters_str = (
            't_0 u_0 t_0_1 u_0_1 t_0_2 u_0_2 t_E t_eff rho rho_1 rho_2 ' +
            't_star t_star_1 t_star_2 pi_E_N pi_E_E s q alpha ds_dt ' +
            'dalpha_dt x_caustic_in x_caustic_out t_caustic_in ' +
            't_caustic_out xi_period xi_semimajor_axis xi_Omega_node ' +
            'xi_inclination xi_argument_of_latitude_reference ' +
            'xi_eccentricity xi_omega_periapsis q_source')
        self._all_MM_parameters = parameters_str.split()
        self._fixed_only_MM_parameters = ['t_0_par', 't_0_xi']

        self._latex_conversion = dict(
            t_0='t_0', u_0='u_0',
            t_0_1='t_{0,1}', u_0_1='u_{0,1}',
            t_0_2='t_{0,2}', u_0_2='u_{0,2}', t_E='t_{\\rm E}',
            t_eff='t_{\\rm eff}', rho='\\rho', rho_1='\\rho_1',
            rho_2='\\rho_2', t_star='t_{\\star}', t_star_1='t_{\\star,1}',
            t_star_2='t_{\\star,2}', pi_E_N='\\pi_{{\\rm E},N}',
            pi_E_E='\\pi_{{\\rm E},E}', s='s', q='q', alpha='\\alpha',
            ds_dt='ds/dt', dalpha_dt='d\\alpha/dt',
            x_caustic_in='x_{\\rm caustic,in}',
            x_caustic_out='x_{\\rm caustic,out}',
            t_caustic_in='t_{\\rm caustic,in}',
            t_caustic_out='t_{\\rm caustic,out}',
            xi_period='\\xi_P',
            xi_semimajor_axis='\\xi_a',
            xi_Omega_node='\\xi_{\\Omega}',
            xi_inclination='\\xi_i',
            xi_argument_of_latitude_reference='\\xi_u',
            xi_eccentricity='\\xi_e',
            xi_omega_periapsis='\\xi_{\\omega}',
            q_source='q_{\\rm source}',
        )

        self._user_parameters = []
        self._other_parameters = []
        self._latex_conversion_user = dict()
        self._latex_conversion_other = dict()

        self._set_default_user_and_other_parameters()

    def _set_default_user_and_other_parameters(self):
        """
        Method to be sub-classed if user defines their own parameters (microlensing or other).
        If you use different microlensing parameters, then define self._user_parameters and
        self._latex_conversion_user.
        If you add non-microlensing parameters (e.g., source distance) then define self._other_parameters and
        self._latex_conversion_other.
        """
        pass

    def _guess_fitting_method(self):
        """
        guess what is the fitting method based on parameters provided
        """
        method = None
        if self._starting_parameters_input is not None:
            method = "EMCEE"
        if self._prior_limits is not None:
            if method is not None:
                raise ValueError(
                    "Both starting_parameters and prior_limits were defined "
                    "which makes impossible to choose the fitting method. "
                    "These settings indicate EMCEE and pyMultiNest "
                    "respectively, and cannot be both set.")
            method = self._guess_MultiNest_or_UltraNest()
        if method is None:
            raise ValueError(
                "No fitting method chosen. You can chose either 'EMCEE' or "
                "'pyMultiNest' and you do it by providing "
                "starting_parameters or prior_limits, respectively.")
        self._fit_method = method

    def _guess_MultiNest_or_UltraNest(self):
        """
        Guess fit_method between MultiNest or UltraNest, based on the
        provided fitting_parameters.
        """
        args_MultiNest = ['basename', 'multimodal', 'evidence tolerance',
                          'n_live_points']
        if all([key in args_MultiNest for key in self._fitting_parameters]):
            return "MultiNest"

        args_UltraNest = ['log directory', 'derived parameter names',
                          'show_status', 'dlogz', 'frac_remain',
                          'max_num_improvement_loops', 'n_live_points']
        if all([key in args_UltraNest for key in self._fitting_parameters]):
            return "UltraNest"

        raise ValueError(
            "Cannot guess fitting method. Provide more parameters in "
            "fitting_parameters.")

    def _check_fitting_method(self):
        """
        Check if fitting method is consistent with the settings.
        """
        if self._fit_method == "emcee":
            self._fit_method = "EMCEE"
            if self._starting_parameters_input is None:
                raise ValueError(
                    "EMCEE fitting method requires starting_parameters.")
        elif self._fit_method in ["multinest", "ultranest"]:
            self._fit_method = self._fit_method.capitalize()
            self._fit_method = self._fit_method.replace("nest", "Nest")
            if self._prior_limits is None:
                msg = "{:} fitting method requires prior_limits."
                raise ValueError(msg.format(self._fit_method))
        else:
            raise ValueError("Invalid fitting method was inserted.")

    def _check_starting_parameters_type(self):
        """
        Check if starting parameters are read from file or
        will be drawn from distributions specified.
        """
        if self._fit_method in ["MultiNest", "UltraNest"]:
            return

        if 'file' in self._starting_parameters_input:
            in_type = 'file'
            keys_expected = {'file', 'parameters'}
            keys = set(self._starting_parameters_input.keys())
            if keys != keys_expected:
                error = ('Wrong starting parameters keys. Expected: ' +
                         str(keys_expected) + '; Provided: ' + str(keys))
                raise KeyError(error)
        else:
            in_type = 'draw'

        self._starting_parameters_type = in_type

    def _set_fit_parameters_unsorted(self):
        """
        Find what are the fitted parameters. It will be sorted later.
        """
        if self._fit_method == "EMCEE":
            if self._starting_parameters_type == 'draw':
                unsorted_keys = self._starting_parameters_input.keys()
            elif self._starting_parameters_type == 'file':
                unsorted_keys = self._get_unsorted_starting_parameters()
            else:
                raise ValueError(
                    'unexpected: ' + str(self._starting_parameters_type))
        elif self._fit_method in ["MultiNest", "UltraNest"]:
            unsorted_keys = self._prior_limits.keys()
        else:
            raise ValueError('unexpected method error')

        self._fit_parameters_unsorted = list(unsorted_keys)
        self._n_fit_parameters = len(self._fit_parameters_unsorted)

    def _get_unsorted_starting_parameters(self):
        """
        Make sure that a list of parameters is provided
        """
        parameters = self._starting_parameters_input['parameters']
        if isinstance(parameters, (str)):
            parameters = parameters.split()

        return parameters

    def _check_imports(self):
        """
        check if all the required packages are imported
        """
        required_packages = set()

        if self._task == 'fit':
            if self._fit_method == 'EMCEE':
                required_packages.add('emcee')
            elif self._fit_method == "MultiNest":
                required_packages.add('pymultinest')
            elif self._fit_method == "UltraNest":
                required_packages.add('ultranest')

            if self._plots is not None and 'triangle' in self._plots:
                required_packages.add('corner')

        if self._plots['best model'] and 'interactive' in self._plots['best model']:
            required_packages.add('plotly')

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
            if "plotly" in failed:
                message += ("\nThe plotly package is required for creating interactive best model plots.")

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
        self._check_other_fit_parameters()
        self._parse_other_output_parameters()
        self._get_datasets()
        self._check_ulens_model_parameters()
        self._get_parameters_ordered()
        self._get_parameters_latex()
        self._set_prior_limits()
        self._parse_fit_constraints()
        if self._fit_method == "EMCEE":
            self._parse_starting_parameters()

        self._check_fixed_parameters()
        self._make_model_and_event()
        self._parse_fitting_parameters()
        if self._fit_method == "EMCEE":
            self._get_starting_parameters()

        self._setup_fit()
        self._run_fit()
        self._finish_fit()
        self._parse_results()
        self._write_residuals()
        self._make_plots()

    def plot_best_model(self):
        """
        Plot the best model.

        The parameters names and their values are taken from __init__()
        keyword ``model``, which is a *dict* and has this information in
        ``model['parameters']`` and ``model['values']``.
        """
        if self._task != "plot":
            raise ValueError('wrong settings to run .plot_best_model()')

        self._check_plots_parameters()
        self._check_model_parameters()
        self._get_datasets()
        self._check_ulens_model_parameters()
        self._check_fixed_parameters()
        self._make_model_and_event()
        self._write_residuals()
        self._make_plots()

    def _check_plots_parameters(self):
        """
        Check if parameters of plots make sense
        """
        allowed_keys = set(['best model', 'trajectory', 'triangle', 'trace'])

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

        if 'trajectory' in self._plots:
            self._check_plots_parameters_trajectory()

        if 'triangle' in self._plots:
            self._check_plots_parameters_triangle()

        if 'trace' in self._plots:
            self._check_plots_parameters_trace()

        names = {key: value.get('file', None)
                 for (key, value) in self._plots.items()}
        done = {}
        for (plot_type, name) in names.items():
            if name is None:
                continue
            if name in done:
                raise ValueError(
                    "Names of output plot files cannot repeat. They repeat "
                    "for: {:} and {:}".format(done[name], plot_type))
            done[name] = plot_type

    def _check_plots_parameters_best_model(self):
        """
        Check if parameters of best model make sense
        """
        allowed = set(['file', 'time range', 'magnitude range', 'legend',
                       'rcParams', 'second Y scale', 'interactive'])
        unknown = set(self._plots['best model'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "best model": {:}'.format(unknown))

        self._set_time_range_for_plot('best model')

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
                    msg = ('The value of {:} (in best model settings)'
                           'must be a dictionary, but you provided {:}')
                    args = [key, type(self._plots['best model'][key])]
                    raise TypeError(msg.format(*args))

        if 'interactive' in self._plots['best model']:
            self._check_plots_parameters_best_model_interactive()

        if 'second Y scale' in self._plots['best model']:
            self._check_plots_parameters_best_model_Y_scale()

    def _set_time_range_for_plot(self, plot_type):
        """
        set time range for best model or triangle plots
        """
        if 'time range' not in self._plots[plot_type]:
            return

        text = self._plots[plot_type]['time range'].split()
        if len(text) != 2:
            msg = ("'time range' for {:} plot should specify 2 values "
                   "(begin and end); got: {:}")
            raise ValueError(
                msg.format(plot_type, self._plots[plot_type]['time range']))
        t_0 = float(text[0])
        t_1 = float(text[1])
        if t_1 < t_0:
            raise ValueError("Incorrect 'time range' for " + plot_type +
                             "plot:\n" + text[0] + " " + text[1])
        self._plots[plot_type]['time range'] = [t_0, t_1]

    def _check_plots_parameters_best_model_interactive(self):
        """
        Check if there is no problem with interactive best plot
        """
        if "second Y scale" in self._plots['best model']:
            msg = "Interactive plot will not have a second Y scale. This feature is not yet implemented."
            raise NotImplementedError(msg)

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
            if settings['magnifications'] != 'optimal':
                raise TypeError(
                    '"best model" -> "second Y scale" -> "magnifications" has '
                    'to be a list or "optimal", not ' +
                    str(type(settings['magnifications'])))
        else:
            for value in settings['magnifications']:
                if not isinstance(value, (int, float)):
                    raise TypeError(
                        'Wrong value in magnifications: ' + str(value))
        if 'labels' not in settings:
            if settings['magnifications'] != 'optimal':
                settings['labels'] = [
                    str(x) for x in settings['magnifications']]
            else:
                settings['labels'] = []
        else:
            if settings['magnifications'] == 'optimal':
                raise ValueError(
                    'In "best model" -> "second Y scale", labels can not be '
                    'provided if "magnifications" is defined as "optimal"')
            if not isinstance(settings['labels'], list):
                raise TypeError(
                    '"best model" -> "second Y scale" -> "labels" has to be '
                    'a list, not ' + str(type(settings['labels'])))
            if len(settings['labels']) != len(settings['magnifications']):
                raise ValueError(
                    'In "best model" -> "second Y scale", labels and '
                    'magnifications must be lists of the same length')

    def _check_plots_parameters_trajectory(self):
        """
        Check if parameters of trajectory plot make sense
        """
        allowed = set(['file', 'time range'])
        unknown = set(self._plots['trajectory'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "trajectory" plot: {:}'.format(unknown))

        self._set_time_range_for_plot('trajectory')

    def _check_plots_parameters_triangle(self):
        """
        Check if parameters of triangle plot make sense
        """
        allowed = set(['file', 'shift t_0'])
        unknown = set(self._plots['triangle'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "triangle" plot: {:}'.format(unknown))

        self._parse_plots_parameter_shift_t_0(self._plots['triangle'])

    def _parse_plots_parameter_shift_t_0(self, settings):
        """
        Check if 'shift t_0' is provided and parse it
        """
        if 'shift t_0' not in settings:
            return

        value = settings['shift t_0']
        if not isinstance(value, bool):
            raise TypeError(
                'For triangle and trace plots, the value of "shift t_0" key '
                'must be of bool type; you provided: ' + str(type(value)))

        self._shift_t_0 = value

    def _check_plots_parameters_trace(self):
        """
        Check if parameters of trace plot make sense
        """
        allowed = set(['file', 'shift t_0'])
        unknown = set(self._plots['trace'].keys()) - allowed
        if len(unknown) > 0:
            raise ValueError(
                'Unknown settings for "trace" plot: {:}'.format(unknown))

        if self._fit_method in ["MultiNest", "UltraNest"]:
            raise ValueError(
                f'Trace plot cannot be requested for {self._fit_method}.')

        self._parse_plots_parameter_shift_t_0(self._plots['trace'])

    def _check_model_parameters(self):
        """
        Check parameters of the MulensModel.Model provided by the user
        directly.
        """
        if self._model_parameters is None:
            self._model_parameters = dict()

        allowed = {
            'coords', 'default method', 'methods', 'methods parameters',
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

        all_parameters = []
        if self._fixed_parameters is not None:
            all_parameters += list(self._fixed_parameters.keys())
        if 'parameters' in keys:
            all_parameters += self._model_parameters['parameters']
        if self._task == 'fit':
            all_parameters += self._fit_parameters_unsorted
        if 'pi_E_E' in all_parameters or 'pi_E_N' in all_parameters:
            if 'coords' not in self._model_parameters:
                raise ValueError("Parallax model requires model['coords'].")

    def _check_other_fit_parameters(self):
        """
        Check if there aren't any other inconsistencies between settings
        """
        if self._fit_method in ["MultiNest", "UltraNest"]:
            if self._min_values is not None or self._max_values is not None:
                msg = "In {:} fitting you cannot set min_values or max_values"
                raise ValueError(msg.format(self._fit_method))

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
                self._parse_other_output_parameters_models(value)
            elif key == 'yaml output':
                self._parse_other_output_parameters_yaml_output(value)
            elif key == 'residuals':
                self._parse_other_output_parameters_residuals(value)
            else:
                raise ValueError('Unrecognized key: ' + str(key) + "\n " +
                                 "Expected keys: models")

    def _parse_other_output_parameters_models(self, values):
        """
        parse information on "other output" -> "models"
        """
        if not isinstance(values, dict):
            raise ValueError('"models" value should also be *dict*, '
                             'got ' + str(type(values)))
        for (key, value) in values.items():
            if key == 'file name':
                self._print_model = True
                self._print_model_i = 0
                if value == '-':
                    self._print_model_file = sys.stdout
                else:
                    try:
                        self._print_model_file = open(value, 'w')
                    except Exception:
                        raise ValueError(
                            'Error while opening file ' + str(value))
            else:
                raise KeyError("Unrecognized key: " + str(key) +
                               "\nExpected keys: 'file name'.")

    def _parse_other_output_parameters_yaml_output(self, values):
        """
        parse information on "other output" -> "yaml output"
        """
        if not isinstance(values, dict):
            raise ValueError('"yaml output" value should also be *dict*, '
                             'got ' + str(type(values)))
        for (key, value) in values.items():
            if key == 'file name':
                self._yaml_results = True
                try:
                    self._yaml_results_file = open(value, 'w')
                except Exception:
                    raise ValueError('Error while opening output '
                                     'YAML file ' + str(value))
                self._yaml_kwargs = {'file': self._yaml_results_file,
                                     'flush': True}
            else:
                raise KeyError("Unrecognized key: " + str(key) +
                               "\nExpected keys: 'file name'.")

    def _parse_other_output_parameters_residuals(self, values):
        """
        parse information on "other output" -> "residuals"
        """
        if not isinstance(values, dict):
            raise ValueError('"residuals" value should also be *dict*, '
                             'got ' + str(type(values)))
        for (key, value) in values.items():
            if key == 'files':
                self._residuals_output = True
                if isinstance(value, list):
                    self._residuals_files = value
                else:
                    if value.count(" ") == 0:
                        self._residuals_files = [value]
                    else:
                        self._residuals_files = value.split()
                self._check_residual_files()
            else:
                raise KeyError("Unrecognized key: " + str(key) +
                               "\nExpected keys: 'file name'.")

    def _check_residual_files(self):
        """
        Check if provided names of output files with residuals make sense.
        We do not check here if the number of files provided is the same
        as the number of input datasets.
        """
        existing = []
        names = []
        for file_name in self._residuals_files:
            if file_name == '-':
                continue
            names.append(file_name)
            if path.exists(file_name):
                existing.append(file_name)
        if len(existing) > 0:
            raise ValueError(
                "Residuals cannot be written to existing files: " +
                str(existing))
        if len(names) != len(set(names)):
            duplicates = set([x for x in names if names.count(x) > 1])
            raise ValueError('some of provided names of residuals files '
                             'repeat: ' + str(duplicates))

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
            dataset = self._get_1_dataset(file_, kwargs)
            self._datasets.append(dataset)

        if self._residuals_output:
            if len(self._residuals_files) != len(self._datasets):
                out = '{:} vs {:}'.format(
                    len(self._datasets), len(self._residuals_files))
                raise ValueError('The number of datasets and files for '
                                 'residuals output do not match: ' + out)

    def _get_1_dataset(self, file_, kwargs):
        """
        Construct a single dataset and possibly rescale uncertainties in it.
        """
        scaling = file_.pop("scale_errorbars", None)
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

        if scaling is not None:
            dataset.scale_errorbars(**scaling)

        return dataset

    def _check_ulens_model_parameters(self):
        """
        Check if there aren't too many parameters.
        Standard check (e.g. if t_0 is defined) are done in mm.Model().
        """
        if self._task == 'fit':
            to_be_checked = set(self._fit_parameters_unsorted)
        elif self._task == 'plot':
            to_be_checked = set(self._model_parameters['parameters'])
        else:
            raise ValueError('unexpected error: ' + str(self._task))
        allowed = self._all_MM_parameters + self._other_parameters + self._user_parameters
        unknown = to_be_checked - set(allowed)
        if len(unknown) > 0:
            raise ValueError('Unknown parameters: {:}'.format(unknown))

    def _get_parameters_ordered(self):
        """
        Order input parameters in some logical way.
        This is useful to make sure the order of printed parameters
        is always the same.
        """
        order = self._all_MM_parameters + self._user_parameters + self._other_parameters
        indexes = sorted(
            [order.index(p) for p in self._fit_parameters_unsorted])

        self._fit_parameters = [order[i] for i in indexes]
        self._fit_parameters_other = [
            order[i] for i in indexes if order[i] in self._other_parameters]
        self._other_parameters_dict = dict()

        if len(self._fit_parameters_other) > 0:
            self._flat_priors = False

    def _get_parameters_latex(self):
        """
        change self._fit_parameters into latex parameters
        """
        conversion = {**self._latex_conversion, **self._latex_conversion_other, **self._latex_conversion_user}

        if self._shift_t_0:
            for key in ['t_0', 't_0_1', 't_0_2']:
                conversion[key] = '\\Delta ' + conversion[key]

        if self._fit_constraints is not None:
            if 'posterior parsing' in self._fit_constraints:
                settings = self._fit_constraints['posterior parsing']
                if 'abs' in settings:
                    if not isinstance(settings['abs'], list):
                        raise ValueError("Error: fit_constraints -> posterior"
                                         " parsing -> abs - list expected")
                    for key in settings['abs']:
                        conversion[key] = "|" + conversion[key] + "|"

        self._fit_parameters_latex = [
            ('$' + conversion[key] + '$') for key in self._fit_parameters]

    def _parse_fitting_parameters(self):
        """
        run some checks on self._fitting_parameters to make sure that
        the fit can be run
        """
        if self._fit_method == 'EMCEE':
            self._parse_fitting_parameters_EMCEE()
            self._get_n_walkers()
        elif self._fit_method == 'MultiNest':
            self._parse_fitting_parameters_MultiNest()
        elif self._fit_method == 'UltraNest':
            self._parse_fitting_parameters_UltraNest()
        else:
            raise ValueError('internal inconsistency')

    def _parse_fitting_parameters_EMCEE(self):
        """
        make sure EMCEE fitting parameters are properly defined
        """
        settings = self._fitting_parameters

        ints_required = ['n_steps']
        required = ints_required

        bools = ['progress']
        ints = ['n_walkers', 'n_burn']
        strings = ['posterior file', 'posterior file fluxes']
        allowed = bools + ints + strings

        self._check_required_and_allowed_parameters(required, allowed)
        self._check_parameters_types(settings, bools=bools,
                                     ints=ints+ints_required, strings=strings)

        self._kwargs_EMCEE = {'initial_state': None,  # It will be set later.
                              'nsteps': self._fitting_parameters['n_steps'],
                              'progress': False}

        if 'progress' in settings:
            self._kwargs_EMCEE['progress'] = settings['progress']

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
            if name[-4:] != '.npy':
                raise ValueError('"posterior file" must end in ".npy", ' +
                                 'got: ' + name)
            if path.exists(name):
                if path.isfile(name):
                    msg = "Existing file " + name + " will be overwritten"
                    warnings.warn(msg)
                else:
                    raise ValueError("The path provided for posterior (" +
                                     name + ") exists and is a directory")
            self._posterior_file_name = name[:-4]
            self._posterior_file_fluxes = None

        if 'posterior file fluxes' in settings:
            fluxes_allowed = ['all', None]
            if settings['posterior file fluxes'] not in fluxes_allowed:
                raise ValueError('Unrecognized "posterior file fluxes": ' +
                                 settings['posterior file fluxes'])
            self._posterior_file_fluxes = settings['posterior file fluxes']

    def _check_required_and_allowed_parameters(self, required, allowed):
        """
        Check if required parameters are provided and there aren't parameters
        that shouldn't be defined.
        """
        settings = self._fitting_parameters
        if settings is None:
            settings = dict()
        full = required + allowed

        for required_ in required:
            if required_ not in settings:
                msg = '{:} method requires fitting parameter: {:}'
                raise ValueError(msg.format(self._fit_method, required_))

        if len(set(settings.keys()) - set(full)) > 0:
            raise ValueError('Unexpected fitting parameters: ' +
                             str(set(settings.keys()) - set(full)))

    def _check_parameters_types(self, settings, bools=None,
                                ints=None, floats=None, strings=None):
        """
        Check if the settings have right type.
        For floats we accept ints as well.
        """
        if bools is None:
            bools = []
        if ints is None:
            ints = []
        if floats is None:
            floats = []
        if strings is None:
            strings = []

        fmt = "For key {:} the expected type is {:}, but got {:}"
        for (key, value) in settings.items():
            if key in bools:
                if not isinstance(value, bool):
                    raise TypeError(
                        fmt.format(key, "bool", str(type(value))))
            elif key in ints:
                if not isinstance(value, int):
                    raise TypeError(
                        fmt.format(key, "int", str(type(value))))
            elif key in floats:
                if not isinstance(value, (float, int)):
                    raise TypeError(
                        fmt.format(key, "float", str(type(value))))
            elif key in strings:
                if not isinstance(value, str):
                    raise TypeError(
                        fmt.format(key, "string", str(type(value))))
            else:
                raise ValueError(
                    "internal bug - no type for key " + key + " specified")

    def _get_n_walkers(self):
        """
        Guess how many walkers (and hence starting values) there will be.
        EMCEE fitting only.
        """
        if self._fit_method != 'EMCEE':
            raise ValueError('internal bug')

        if 'n_walkers' in self._fitting_parameters:
            self._n_walkers = self._fitting_parameters['n_walkers']
        else:
            if self._starting_parameters_type == 'draw':
                self._n_walkers = 4 * self._n_fit_parameters
            elif self._starting_parameters_type != 'file':
                raise ValueError(
                    'Unexpected: ' + self._starting_parameters_type)

    def _parse_fitting_parameters_MultiNest(self):
        """
        make sure MultiNest fitting parameters are properly defined
        """
        self._kwargs_MultiNest = dict()

        settings = self._fitting_parameters
        if settings is None:
            settings = dict()

        required = []
        bools = ['multimodal']
        ints = ['n_live_points']
        strings = ['basename']
        floats = ['sampling efficiency', 'evidence tolerance']
        allowed = strings + bools + ints + floats

        self._check_required_and_allowed_parameters(required, allowed)
        self._check_parameters_types(settings, bools, ints, floats, strings)

        self._kwargs_MultiNest['multimodal'] = False

        keys = {"basename": "outputfiles_basename",
                "sampling efficiency": "sampling_efficiency",
                "evidence tolerance": "evidence_tolerance"}
        same_keys = ["multimodal", "n_live_points"]
        keys = {**keys, **{key: key for key in same_keys}}

        self._set_dict_safely(self._kwargs_MultiNest, settings, keys)

        self._kwargs_MultiNest['importance_nested_sampling'] = (
            not self._kwargs_MultiNest['multimodal'])
        self._MN_temporary_files = False
        if 'basename' not in settings:
            print("No base for MultiNest output provided.")
            self._kwargs_MultiNest['outputfiles_basename'] = (
                tempfile.mkdtemp('_MM_ex16_pyMN') + sep)
            self._MN_temporary_files = True
        self._check_output_files_MultiNest()

    def _set_dict_safely(self, target, source, keys_mapping):
        """
        For each key in keys_mapping (*dict*) check if it is in
        source (*dict*). If it is, then set
        target[keys_mapping[key]] to source[key].
        """
        for (key_in, key_out) in keys_mapping.items():
            if key_in in source:
                target[key_out] = source[key_in]

    def _check_output_files_MultiNest(self):
        """
        Check if output files exist and warn about overwriting them.

        If they directory doesn't exist then raise error.
        """
        root = self._kwargs_MultiNest['outputfiles_basename']
        if len(path.dirname(root)) != 0 and not path.isdir(path.dirname(root)):
            msg = 'directory for output files does not exist; root path: '
            raise ValueError(msg + root)
        if path.isdir(root):
            root += sep

        check = (
            'resume.dat phys_live.points live.points phys_live-birth.txt '
            'ev.dat dead-birth.txt .txt stats.dat post_equal_weights.dat'
            'post_separate_strict.dat post_separate.dat summary.txt').split()
        if self._kwargs_MultiNest['importance_nested_sampling']:
            check += ['IS.points', 'IS.ptprob', 'IS.iterinfo']

        root = self._kwargs_MultiNest['outputfiles_basename']
        if path.isdir(root):
            root += sep

        existing = []
        for check_ in check:
            file_name = root + check_
            if path.isfile(file_name):
                existing.append(file_name)

        if len(existing) > 0:
            message = "\n\n Existing files will be overwritten "
            message += "(unless you kill this process)!!!\n"
            warnings.warn(message + str(existing) + "\n")

    def _parse_fitting_parameters_UltraNest(self):
        """
        Make sure UltraNest fitting parameters are properly defined
        """
        self._kwargs_UltraNest = dict()
        self._kwargs_UltraNest['viz_callback'] = False

        settings = self._fitting_parameters.copy()
        if settings is None:
            settings = dict()

        required = []
        bools = ['show_status']
        ints = ['min_num_live_points', 'max_num_improvement_loops']
        if 'n_live_points' in settings:
            ints[0] = 'n_live_points'
        strings = ['log directory', 'derived parameter names']
        floats = ['dlogz', 'frac_remain']
        allowed = strings + bools + ints + floats

        self._check_required_and_allowed_parameters(required, allowed)
        self._check_parameters_types(settings, bools, ints, floats, strings)
        self._log_dir_UltraNest = settings.pop("log directory", None)
        value = settings.pop("derived parameter names", "")
        self._derived_param_names_UltraNest = value.split()
        self._check_dir_and_parameter_names_Ultranest()

        keys = {"n_live_points": "min_num_live_points"}
        same_keys = ["min_num_live_points", 'max_num_improvement_loops', "show_status", "dlogz", "frac_remain"]
        keys = {**keys, **{key: key for key in same_keys}}
        self._set_dict_safely(self._kwargs_UltraNest, settings, keys)

    def _check_dir_and_parameter_names_Ultranest(self):
        """
        Checks if the path to `log directory` exists and is a directory,
        and also if the number of `derived parameter names` matches the
        number of derived fluxes.
        """
        if self._log_dir_UltraNest is not None:
            if not path.exists(self._log_dir_UltraNest):
                raise ValueError("log directory value in fitting_parameters does not exist.")
            elif not path.isdir(self._log_dir_UltraNest):
                raise ValueError("log directory value in fitting_parameters exists, but it is a file.")

        n_datasets = len(self._datasets)
        n_fluxes = self._n_fluxes_per_dataset
        if len(self._derived_param_names_UltraNest) != n_datasets * n_fluxes:
            raise ValueError("The number of `derived parameter names` must match the number of derived fluxes.")

    def _set_prior_limits(self):
        """
        Set minimum and maximum values of the prior space
        """
        if self._fit_method == 'EMCEE':
            self._set_prior_limits_EMCEE()
        elif self._fit_method in ['MultiNest', 'UltraNest']:
            self._set_prior_limits_MultiNest()
        else:
            raise ValueError('internal bug')

    def _set_prior_limits_EMCEE(self):
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
                    fmt = ("This doesn't make sense: for {:} the lower limit "
                           "is larger than the upper limit: {:} vs {:}")
                    raise ValueError(fmt.format(key, self._min_values[key], self._max_values[key]))

        self._min_values_indexed = self._parse_min_max_values_single(self._min_values)
        self._max_values_indexed = self._parse_min_max_values_single(self._max_values)

    def _parse_min_max_values_single(self, limits):
        """
        change dict that has str as key to index as key
        """
        out = dict()
        if len(limits) == 0:
            return out

        for (key, value) in limits.items():
            if key not in self._fit_parameters:
                fmt = 'Key provided in limits: {:}\nis not one of the parameters for fitting: {:}'
                raise ValueError(fmt.format(key, self._fit_parameters))

            index = self._fit_parameters.index(key)
            out[index] = value

        return out

    def _set_prior_limits_MultiNest(self):
        """
        Set prior limits and transformation constants (from unit cube to MM parameters).
        """
        min_values = []
        max_values = []
        for parameter in self._fit_parameters:
            if parameter not in self._prior_limits:
                raise ValueError("internal issue")
            values = self._prior_limits[parameter]
            if isinstance(values, str):
                values = values.split()
            if not isinstance(values, list) or len(values) != 2:
                raise ValueError(
                    "prior_limits for " + parameter + " could not be "
                    "processed: " + str(self._prior_limits[parameter]))
            try:
                values = [float(v) for v in values]
            except Exception:
                raise ValueError(
                    "couldn't get floats for prior_limits of " + parameter +
                    ": " + str(self._prior_limits[parameter]))
            if values[0] >= values[1]:
                raise ValueError(
                    "This won't work - wrong order of limits for " +
                    parameter + ": " + str(values))
            min_values.append(values[0])
            max_values.append(values[1])

        self._min_values = np.array(min_values)
        self._max_values = np.array(max_values)
        self._range_values = self._max_values - self._min_values

    def _parse_fit_constraints(self):
        """
        Parse the fitting constraints that are not simple limits on parameters
        """
        self._prior_t_E = None
        self._priors = None

        if self._fit_constraints is None:
            self._set_default_fit_constraints()
            return

        self._check_fit_constraints()
        self._parse_fit_constraints_keys()
        self._parse_fit_constraints_fluxes()
        self._parse_fit_constraints_posterior()

        if 'prior' in self._fit_constraints:
            self._parse_fit_constraints_prior()

    def _check_fit_constraints(self):
        """
        Run checks on self._fit_constraints
        """
        if self._fit_method == 'MultiNest':
            raise NotImplementedError(
                "Currently no fit_constraints are implemented for MultiNest "
                "fit. Please contact Radek Poleski with a specific request.")

        if self._fit_method == 'UltraNest':
            allowed_keys = {'negative_blending_flux_sigma_mag', 'prior'}
            used_keys = set(self._fit_constraints.keys())
            if not used_keys.issubset(allowed_keys):
                raise NotImplementedError(
                    "The supported fit_constraints options for UltraNest are"
                    " `negative_blending_flux_sigma_mag` and `prior`.")
            if 'prior' in used_keys:
                if set(self._fit_constraints['prior'].keys()) != {'t_E'}:
                    raise ValueError(
                        "Only `t_E` is allowed in fit_constraints['prior'].")

        if isinstance(self._fit_constraints, list):
            raise TypeError(
                "In version 0.5.0 we've changed type of 'fit_constraints' " +
                "from list to dict. Please correct your input and re-run " +
                "the code. Most probably what you need is:\n" +
                "fit_constraints = {'no_negative_blending_flux': True}")

    def _parse_fit_constraints_keys(self):
        """
        Validate the keys in the provided fit_constraints.
        """
        allowed_keys_flux = {
            "no_negative_blending_flux", "negative_blending_flux_sigma_mag"}
        allowed_keys_color = {'color', 'color source 1', 'color source 2'}
        allowed_keys = {*allowed_keys_flux, *allowed_keys_color,
                        "prior", "posterior parsing"}

        used_keys = set(self._fit_constraints.keys())
        if len(used_keys - allowed_keys) > 0:
            raise ValueError('unrecognized constraint: {:}'.format(
                used_keys - allowed_keys))
        if len(used_keys.intersection(allowed_keys_flux)) == 2:
            raise ValueError(
                'you cannot specify both no_negative_blending_flux and ' +
                'negative_blending_flux_sigma_mag')
        if "no_negative_blending_flux" not in self._fit_constraints:
            self._fit_constraints["no_negative_blending_flux"] = False

        self._check_color_constraints_conflict(allowed_keys_color)

    def _set_default_fit_constraints(self):
        """
        Set default fitting constraints if none are provided.
        """
        self._fit_constraints = {"no_negative_blending_flux": False}
        self._parse_posterior_abs = list()

    def _check_color_constraints_conflict(self, allowed_keys_color):
        """
        Check for conflicts among color constraints.
        """
        used_keys = set(self._fit_constraints.keys())

        if len(used_keys.intersection(allowed_keys_color)) >= 2:
            if 'color' in used_keys:
                raise ValueError(
                    'You cannot specify both color and ' +
                    str(used_keys.intersection(allowed_keys_color)-{'color'}))

    def _parse_fit_constraints_fluxes(self):
        """
        Process each constraint fit_constraints.
        """
        for key, value in self._fit_constraints.items():
            if key == "negative_blending_flux_sigma_mag":
                self._parse_fit_constraints_soft_blending(key, value)
            elif key in ['color', 'color source 1', 'color source 2']:
                self._parse_fit_constraints_color(key, value)

    def _parse_fit_constraints_soft_blending(self, key, value):
        """
        Check if soft fit constraint on blending flux are correctly defined.
        """
        if isinstance(value,  float):
            sigma = float(value)
            sets = list(range(len(self._datasets)))

        else:

            sigma = float(value.split()[0])
            sets = list(map(self._get_no_of_dataset,
                            shlex.split(value, posix=False)[1:]))
            if len(sets) > len(self._datasets):
                raise ValueError(
                    "dataset number specified in" +
                    "negative_blending_flux_sigma_mag" +
                    "do not match with provided datasets")

        self._fit_constraints[key] = [
            mm.Utils.get_flux_from_mag(sigma), sets]

    def _parse_fit_constraints_color(self, key, value):
        """
        Check if fit constraint on color are correctly defined.
        """
        self._check_unique_datasets_labels()
        words = shlex.split(value, posix=False)

        if len(words) != 5 or words[0] != 'gauss':
            msg = "Something went wrong in parsing prior for "
            msg += "{:}: {:}"
            if len(words) == 3 and words[0] == 'gauss':
                msg += "color priors require the specification"
                msg += "of datasets that should be used for color calculation"
            raise ValueError(msg.format(key, value))
        try:
            settings = [words[0], float(words[1]), float(
                words[2]), words[3], words[4]]
        except Exception:
            raise ValueError('error in parsing: ' + words[1] + " " +
                             words[2] + " " + words[3] + " " + words[4])
        if settings[2] < 0.:
            raise ValueError('sigma cannot be negative: ' + words[2])

        settings[3] = self._get_no_of_dataset(settings[3])
        settings[4] = self._get_no_of_dataset(settings[4])

        if settings[3] == settings[4]:
            raise ValueError(
                "in " + key + " color have to be from different datasets")
        n = len(self._datasets)-1
        if (0 >= settings[3] >= n) or (0 >= settings[4] >= n):
            raise ValueError(
                "label specified in color prior" +
                "do not match with provided datasets")

        self._fit_constraints[key] = settings

    def _check_unique_datasets_labels(self):
        """
        Check if the labels of datasets are unique.
        """
        labels = [
            dataset.plot_properties['label']
            for dataset in self._datasets
        ]
        if len(labels) != len(set(labels)):
            raise ValueError("Declared labels of datasets must be unique.")

    def _parse_fit_constraints_prior(self):
        """
        Check if priors in fit constraint are correctly defined.
        """
        priors = dict()
        for (key, value) in self._fit_constraints['prior'].items():
            if key == 't_E':
                if value in ["Mroz et al. 2017", "Mroz et al. 2020"]:
                    self._prior_t_E = value.replace(" et al. 20", "+")
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

    def _parse_fit_constraints_posterior(self):
        """
        Parse constraints on what is done with posterior.
        """
        if 'posterior parsing' not in self._fit_constraints:
            self._parse_posterior_abs = list()
            return

        if self._fit_method != "EMCEE":
            raise ValueError('Input in "posterior parsing" is allowed only'
                             ' for EMCEE')

        allowed_keys = {"abs"}
        settings = self._fit_constraints['posterior parsing']
        unknown = set(settings.keys()) - allowed_keys
        if len(unknown) > 0:
            msg = "Unrecognized key in fit_constraints -> 'posterior parsing':"
            raise KeyError(msg + ' ' + str(unknown))

        if 'abs' in settings:
            self._parse_posterior_abs = settings['abs']
            for parameter in self._parse_posterior_abs:
                if parameter not in self._fit_parameters:
                    raise ValueError(
                        "Error - you can calculate absolute value only of "
                        "a parameter which is fitted, not: " + parameter)

    def _get_no_of_dataset(self, label):
        """
        Returns the index of a dataset with a specific label.
        Parameters :
            label: *str* ,*int*
              Label of the dataset defined by
              MulensData.plot_properties['label'],
              or name of the data file if label is not specified,
              or a sequential index of the dataset.

        Returns :
          index: *int*
          Sequential index of the dataset from [0,1,...,n_datasets-1]

        """

        if '"' in label:
            label = label.strip('"')
        for (i, dataset) in enumerate(self._datasets):
            if dataset.plot_properties['label'] == label:
                return i
        raise KeyError(
            "Unrecognized dataset label in fit_constraints/prior: " + label)

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
        elif self._prior_t_E == 'Mroz+20':
            # XXX - TO DO:
            # - documentation
            # - smooth the input data from M+20 and note that
            x = np.array([
                0.74, 0.88, 1.01, 1.15, 1.28, 1.42, 1.55, 1.69,
                1.82, 1.96, 2.09, 2.23, 2.36, 2.50, 2.63])
            y = np.array([
                82.04, 94.98, 167.76, 507.81, 402.08, 681.61, 1157.51,
                1132.80, 668.12, 412.20, 236.14, 335.34, 74.88, 52.64, 97.78])
            dx = (x[1] - x[0]) / 2.
            x_min = x[0] - dx
            x_max = x[-1] + dx
            function = interp1d(x, np.log(y),
                                kind='cubic', fill_value="extrapolate")
        else:
            raise ValueError('unexpected internal error')

        self._prior_t_E_data['x_min'] = x_min
        self._prior_t_E_data['x_max'] = x_max
        self._prior_t_E_data['y_min'] = function(x_min)
        self._prior_t_E_data['y_max'] = function(x_max)
        self._prior_t_E_data['function'] = function

    def _parse_starting_parameters(self):
        """
        Read the starting parameters from file or
        change the format of provided information.
        """
        if self._starting_parameters_type == 'file':
            self._read_starting_parameters_from_file()
        elif self._starting_parameters_type == 'draw':
            self._parse_starting_parameters_to_be_drawn()
        else:
            raise ValueError(
                'unexpected: ' + str(self._starting_parameters_type))

    def _read_starting_parameters_from_file(self):
        """
        Read starting parameters from a file.
        """
        file_name = self._starting_parameters_input['file']
        try:
            data = np.loadtxt(file_name, unpack=True)
        except Exception:
            raise RuntimeError('Error while reading file: ' + file_name)

        if len(data.shape) != 2 or data.shape[0] != self._n_fit_parameters:
            raise ValueError(
                'Something wrong in shape of data read from ' + file_name)

        self._starting_parameters_list = data.T.tolist()
        self._n_walkers = len(self._starting_parameters_list)

    def _parse_starting_parameters_to_be_drawn(self):
        """
        replace self._starting_parameters with dict that has values
        [*str*, *float*, ...]
        and make basic checks
        """
        accepted_types = ['gauss', 'uniform', 'log-uniform']

        out = dict()
        for (key, value) in self._starting_parameters_input.items():
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

        fixed = set(self._fixed_parameters.keys())

        allowed = set(self._all_MM_parameters +
                      self._fixed_only_MM_parameters +
                      self._other_parameters)
        unknown = fixed - allowed
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
        for key in self._other_parameters:
            try:
                parameters.pop(key)
            except Exception:
                pass

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
                model.default_magnification_method = \
                    self._model_parameters['default method']
            if 'methods' in self._model_parameters:
                model.set_magnification_methods(
                    self._model_parameters['methods'])
            if 'methods parameters' in self._model_parameters:
                model.set_magnification_methods_parameters(
                    self._model_parameters['methods parameters'])
            if 'methods source 1' in self._model_parameters:
                self._model.set_magnification_methods(
                    self._model_parameters['methods source 1'], 1)
            if 'methods source 2' in self._model_parameters:
                self._model.set_magnification_methods(
                    self._model_parameters['methods source 2'], 2)

        self._event = mm.Event(self._datasets, self._model)
        self._event.sum_function = 'numpy.sum'

        self._set_n_fluxes()

    def _set_n_fluxes(self):
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
        if self._fixed_parameters is None:
            parameters = dict()
        else:
            parameters = {**self._fixed_parameters}

        if self._task == 'plot':
            parameters.update(dict(zip(
                self._model_parameters['parameters'],
                self._model_parameters['values'])))
            # XXX this is some kind of a hack:
            self._best_model_theta = []
            self._fit_parameters = []
        elif self._task == 'fit':
            if self._fit_method == 'EMCEE':
                parameters.update(self._get_example_parameters_EMCEE())
            elif self._fit_method in ['MultiNest', 'UltraNest']:
                means = 0.5 * (self._max_values + self._min_values)
                parameters.update(dict(zip(self._fit_parameters, means)))
                if "x_caustic_in" in self._fit_parameters:
                    index = self._fit_parameters.index("x_caustic_in")
                    parameters["x_caustic_in"] = (
                        self._min_values[index] +
                        np.random.uniform() * self._range_values[index])
            else:
                raise ValueError('internal value')
        else:
            raise ValueError('internal value')

        parameters = self._transform_parameters(parameters)
        return parameters

    def _transform_parameters(self, parameters):
        """
        Method to be sub-classed if user defines their own microlensing parameters.
        It takes a dict as input and returns a dict.
        """
        return parameters

    def _get_example_parameters_EMCEE(self):
        """
        get sample values of parameters for EMCEE - only to make mm.Model
        """
        if self._starting_parameters_type == 'file':
            return self._extract_example_parameters_EMCEE()
        elif self._starting_parameters_type == 'draw':
            return self._draw_example_parameters_EMCEE()
        else:
            raise ValueError(
                'unexpected: ' + str(self._starting_parameters_type))

    def _extract_example_parameters_EMCEE(self):
        """
        Provide example parameters that will be used to initialize
        MM.Model and MM.Event. Sets of parameters are already read,
        so we just take 0-th set.
        """
        keys = self._fit_parameters_unsorted
        values = self._starting_parameters_list[0]
        return dict(zip(keys, values))

    def _draw_example_parameters_EMCEE(self):
        """
        Randomly draw parameters or EMCEE - only to make mm.Model.
        """
        parameters = dict()
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
                elif value[0] in ['uniform', 'log-uniform']:
                    parameters[key] = (value[1] + value[2]) / 2.
                else:
                    raise ValueError('internal error: ' + value[0])
        return parameters

    def _get_starting_parameters(self):
        """
        Generate random parameters or read them from file.
        Check if there are enough of them and store.
        """
        if self._starting_parameters_type == 'file':
            starting = self._starting_parameters_list
        elif self._starting_parameters_type == 'draw':
            starting = self._generate_random_parameters()
        else:
            raise ValueError(
                'unexpected: ' + str(self._starting_parameters_type))

        self._check_and_store_generated_random_parameters(starting)

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

        return np.array(starting).T.tolist()

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

    def _check_and_store_generated_random_parameters(self, starting):
        """
        Check if the set of points provided has at least self._n_walkers
        points inside the prior and save these.
        """
        verified = self._check_generated_random_parameters(starting)

        if len(verified) < self._n_walkers:
            raise ValueError(
                "Couldn't generate required starting points in a prior. "
                "Most probably you have to correct at least one of: "
                "starting_parameters, min_values, max_values, or "
                "fit_constraints.\nGot " + str(len(verified)) + " walkers in "
                "the prior, but required " + str(self._n_walkers) + ".\n"
                "If you think the code should work with your settings, "
                "then please contact Radek Poleski.")

        self._kwargs_EMCEE['initial_state'] = verified

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

        return out

    def _ln_prob(self, theta):
        """
        Log probability of the model - combines _ln_prior(), _ln_like(),
        and constraints which include fluxes.

        NOTE: we're using np.log(), i.e., natural logarithms.
        """
        if self._fit_method == "EMCEE":
            ln_prior = self._ln_prior(theta)
            if not np.isfinite(ln_prior):
                return self._return_ln_prob(-np.inf)
        elif self._fit_method == "UltraNest":
            ln_prior = self._ln_prior_t_E() if self._prior_t_E else 0.

        ln_like = self._ln_like(theta)
        if not np.isfinite(ln_like):
            return self._return_ln_prob(-np.inf)

        ln_prob = ln_prior + ln_like

        fluxes = self._get_fluxes()

        ln_prior_flux = self._run_flux_checks_ln_prior(fluxes)
        if not np.isfinite(ln_prior_flux):
            return self._return_ln_prob(-np.inf)

        ln_prob += ln_prior_flux

        if self._fit_method == "EMCEE":
            self._update_best_model_EMCEE(ln_prob, theta, fluxes)
            return self._return_ln_prob(ln_prob, fluxes)

        return ln_prob

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
        if self._task == 'plot':
            return

        parameters = dict(zip(self._fit_parameters, theta))
        parameters = self._transform_parameters(parameters)

        if len(self._fit_parameters_other) == 0:
            for (parameter, value) in parameters.items():
                setattr(self._model.parameters, parameter, value)
        else:
            for (parameter, value) in parameters.items():
                if parameter not in self._fit_parameters_other:
                    setattr(self._model.parameters, parameter, value)
                else:
                    self._other_parameters_dict[parameter] = value

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

        if "x_caustic_in" in self._model.parameters.parameters:
            self._set_model_parameters(theta)
            if not self._check_valid_Cassan08_trajectory():
                return outside

        ln_prior = inside

        if self._prior_t_E is not None:
            self._set_model_parameters(theta)
            ln_prior += self._ln_prior_t_E()

        if self._priors is not None:
            self._set_model_parameters(theta)
            for (parameter, prior_settings) in self._priors.items():
                if parameter in ['pi_E_N', 'pi_E_E']:
                    # Other parameters can be added here. XXX
                    value = self._model.parameters.parameters[parameter]
                    ln_prior += self._get_ln_prior_for_1_parameter(
                        value, prior_settings)

                else:
                    raise ValueError('prior not handled: ' + parameter)

        return ln_prior

    def _check_valid_Cassan08_trajectory(self):
        """
        Check if current model (that has to be in Cassan08 parameterization)
        has valid trajectory.
        """
        sampling = self._model.parameters.uniform_caustic_sampling
        valid = sampling.check_valid_trajectory(
            self._model.parameters.x_caustic_in,
            self._model.parameters.x_caustic_out)
        return valid

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
            raise ValueError('unexpected internal error: ' + self._prior_t_E)

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
        out = -0.5 * chi2

        if self._print_model:
            self._print_current_model(theta, chi2)

        if self._task == 'fit' and len(self._other_parameters_dict) > 0:
            out += self._get_ln_probability_for_other_parameters()

        return out

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

    def _get_ln_probability_for_other_parameters(self):
        """
        Function that defines calculation of
        ln(probability(other_parameters)).
        The logarithm is to the base of the constant e.

        If you have defined other_parameters, then you have to implement
        child class of UlensModelFit and re-define this function.

        NOTE: your implementation should primarily use *dict*
        `self._other_parameters_dict`, but all MM parameters are already
        set and kept in *MM.Event* instance `self._event`.
        """
        raise NotImplementedError(
            'You have to subclass this class and re-implement this function')

    def _get_fluxes(self):
        """
        Extract all fluxes and return them in a list.
        """
        fluxes = []
        for (i, dataset) in enumerate(self._datasets):
            fluxes += self._event.fits[i].source_fluxes.tolist()
            fluxes.append(self._event.fits[i].blend_flux)

        return fluxes

    def _sumup_inside_prior(self, fluxes, key, inside, index_plus):
        """
        Calculates the contribution to the ln_prior
        from specified color constraints
        Parameters :
            fluxes: *array*
                Array with fluxes of the current model.
            key: *str*
                constrain key.
            inside: *float*
                ln_prior contribution
            index_plus: *int*
                For a single source, index_plus=0;
                for a binary source, index_plus=0 or 1.
        Returns :
            inside: *float*
                Evaluated ln_prior contribution
        """
        settings = self._fit_constraints[key]
        index1 = (settings[3])*self._n_fluxes_per_dataset + index_plus
        index2 = (settings[4])*self._n_fluxes_per_dataset + index_plus
        value = mm.Utils.get_mag_from_flux(
            fluxes[index1])-mm.Utils.get_mag_from_flux(fluxes[index2])
        inside += self._get_ln_prior_for_1_parameter(value, settings[:-2])

        return inside

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

        inside += self._apply_negative_blending_flux_sigma_mag_prior(fluxes)
        if self._fit_method == "EMCEE":
            inside += self._apply_color_prior(fluxes)

        return inside

    def _apply_negative_blending_flux_sigma_mag_prior(self, fluxes):
        """
        Apply the negative blending flux sigma magnitude priotr.
        """
        inside = 0.0
        key = "negative_blending_flux_sigma_mag"

        if key in self._fit_constraints:
            sigma, datasets = self._fit_constraints[key]
            for i, dataset in enumerate(self._datasets):
                if i in datasets:
                    blend_index = ((i + 1) * self._n_fluxes_per_dataset) - 1
                    if fluxes[blend_index] < 0.0:
                        inside += -0.5 * (fluxes[blend_index] / sigma) ** 2

        return inside

    def _apply_color_prior(self, fluxes):
        """
        Apply the color constraints.
        """
        inside = 0.0
        key = 'color'
        if key in self._fit_constraints:
            for i in range(self._n_fluxes_per_dataset - 1):
                inside += self._sumup_inside_prior(fluxes, key, inside, i)

        key = 'color source 1'
        if key in self._fit_constraints:
            inside += self._sumup_inside_prior(fluxes, key, inside, 0)

        key = 'color source 2'
        if key in self._fit_constraints:
            inside += self._sumup_inside_prior(fluxes, key, inside, 1)
        return inside

    def _update_best_model_EMCEE(self, ln_prob, theta, fluxes):
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
        if self._fit_method == 'EMCEE':
            self._setup_fit_EMCEE()
        elif self._fit_method == 'MultiNest':
            self._setup_fit_MultiNest()
        elif self._fit_method == 'UltraNest':
            self._setup_fit_UltraNest()
        else:
            raise ValueError('internal bug')

    def _setup_fit_EMCEE(self):
        """
        Setup EMCEE fit
        """
        self._sampler = emcee.EnsembleSampler(
            self._n_walkers, self._n_fit_parameters, self._ln_prob)

    def _setup_fit_MultiNest(self):
        """
        Prepare MultiNest fit
        """
        self._kwargs_MultiNest['Prior'] = self._transform_unit_cube
        self._kwargs_MultiNest['LogLikelihood'] = self._ln_like_MN
        self._kwargs_MultiNest['resume'] = False

        self._kwargs_MultiNest['n_dims'] = self._n_fit_parameters
        self._kwargs_MultiNest['n_params'] = self._kwargs_MultiNest['n_dims']
        if self._return_fluxes:
            self._kwargs_MultiNest['n_params'] += self._n_fluxes

    def _setup_fit_UltraNest(self):
        """
        Prepare UltraNest fit, declaring sampler instance.
        If the names of the derived parameters are not given, the source
        and blending fluxes are assigned, with indexes depending on the
        number of sources (s1, s2) and datasets (_1, _2).
        """
        if self._return_fluxes:
            if len(self._derived_param_names_UltraNest) == 0:
                if self._flux_names is None:
                    self._flux_names = self._get_fluxes_names_to_print()
                self._derived_param_names_UltraNest = self._flux_names

        n_dims = self._n_fit_parameters
        n_params = n_dims + self._n_fluxes_per_dataset * self._return_fluxes
        t_kwargs = {'n_dims': n_dims, 'n_params': n_params}
        self._sampler = ultranest.ReactiveNestedSampler(
            self._fit_parameters, self._ln_prob,
            transform=lambda cube: self._transform_unit_cube(cube, **t_kwargs),
            derived_param_names=self._derived_param_names_UltraNest,
            log_dir=self._log_dir_UltraNest
        )

    def _transform_unit_cube(self, cube, n_dims, n_params):
        """
        Transform MultiNest/UltraNest unit cube to microlensing parameters.

        MultiNest: based on SafePrior() in
        https://github.com/JohannesBuchner/PyMultiNest/blob/master/
        pymultinest/solve.py
        UltraNest: based on the above and UltraNest documentation
        https://johannesbuchner.github.io/UltraNest/example-sine-line.html

        NOTE: We call self._ln_like() here (and remember the result)
        because in MultiNest you can add fluxes only in "prior" function,
        not in likelihood function.
        """
        cube_out = self._min_values + cube[:n_dims] * self._range_values
        if self._fit_method == "MultiNest":
            for i in range(n_dims):
                cube[i] = cube_out[i]

        if "x_caustic_in" in self._model.parameters.parameters:
            self._set_model_parameters(cube_out)
            if not self._check_valid_Cassan08_trajectory():
                self._last_ln_like = -1.e300
                self._last_theta = cube_out
                if self._return_fluxes:
                    self._last_fluxes = np.zeros(n_params - n_dims)
                    if self._fit_method == "MultiNest":
                        for i in range(n_dims, n_params):
                            cube[i] = 0.
                        return
                    cube_out = np.append(cube_out, self._last_fluxes)
                    return cube_out

        self._last_ln_like = self._ln_like(cube_out)
        self._last_theta = cube_out

        if self._return_fluxes:
            fluxes = self._get_fluxes()
            self._last_fluxes = fluxes
            if self._fit_method == "UltraNest":
                cube_out = np.append(cube_out, fluxes)
                return cube_out
            for i in range(n_dims, n_params):
                cube[i] = fluxes[i-n_dims]

    def _ln_like_MN(self, theta, n_dim, n_params, lnew):
        """
        Calculate likelihood and save if its best model.
        This is used for MultiNest fitting.
        """
        for i in range(n_dim):
            if self._last_theta[i] != theta[i]:
                msg = "internal bug:\n{:}\n{:}\n{:}"
                raise ValueError(msg.format(i, self._last_theta[i], theta[i]))

        ln_like = self._last_ln_like

        if ln_like > self._best_model_ln_prob:
            self._best_model_ln_prob = ln_like
            self._best_model_theta = np.array(theta[:n_dim])
            if self._return_fluxes:
                self._best_model_fluxes = self._last_fluxes

        ln_max = -1.e300
        if not np.isfinite(ln_like) or ln_like < ln_max:
            if not np.isfinite(ln_like):
                msg = "problematic likelihood: {:}\nfor parameters: {:}"
                warnings.warn(msg.format(ln_like, theta))
            ln_like = ln_max

        return ln_like

    def _run_fit(self):
        """
        Call the method that does the fit.
        """
        if self._fit_method == 'EMCEE':
            self._run_fit_EMCEE()
        elif self._fit_method == 'MultiNest':
            self._run_fit_MultiNest()
        elif self._fit_method == 'UltraNest':
            self._run_fit_UltraNest()
        else:
            raise ValueError('internal bug')

    def _run_fit_EMCEE(self):
        """
        Run EMCEE
        """
        self._sampler.run_mcmc(**self._kwargs_EMCEE)

    def _run_fit_MultiNest(self):
        """
        Run MultiNest fit
        """
        mn_run(**self._kwargs_MultiNest)

    def _run_fit_UltraNest(self):
        """
        Run Ultranest fit
        """
        min_n_live = self._kwargs_UltraNest.get("min_num_live_points", 400)
        cluster_n_live = 40 if min_n_live >= 40 else min_n_live
        self._kwargs_UltraNest['cluster_num_live_points'] = cluster_n_live

        self._result_UltraNest = self._sampler.run(**self._kwargs_UltraNest)

    def _finish_fit(self):
        """
        Make the things that are necessary after the fit is done.
        Currently it's just closing the file with all models and
        reads the output files (only for MultiNest).
        """
        if self._print_model:
            if self._print_model_file is not sys.stdout:
                self._print_model_file.close()
                self._print_model = False

        if self._fit_method == 'MultiNest':
            base = self._kwargs_MultiNest['outputfiles_basename']
            self._analyzer = Analyzer(n_params=self._n_fit_parameters,
                                      outputfiles_basename=base)
            self._analyzer_data = self._analyzer.get_data()

            if self._kwargs_MultiNest['multimodal']:
                self._read_multimode_posterior_MultiNest()
                self._read_multimode_best_models_MultiNest()

    def _read_multimode_posterior_MultiNest(self):
        """
        We read data from MultiNest output file [root]post_separate.dat.
        It has 2 empty lines then info on a mode at this repeats N_modes
        times. We read it twice - first to get all the data and then to
        find how many samples there are in each mode.
        """
        data = np.loadtxt(self._analyzer.post_file)

        n_samples = []
        skip = False
        with open(self._analyzer.post_file) as in_data:
            for line in in_data.readlines():
                if skip:
                    skip = False
                    continue
                if len(line) == 1:
                    skip = True
                    n_samples.append(0)
                else:
                    n_samples[-1] += 1

        if data.shape[0] != sum(n_samples):
            raise ValueError('Error in file ' + self._analyzer.post_file +
                             str(data.shape) + " vs. " + str(n_samples))

        self._MN_samples_modes_all = data
        self._MN_modes_indexes = n_samples
        self._n_modes = len(n_samples)

    def _read_multimode_best_models_MultiNest(self):
        """
        read and process information about best models
        (i.e., highest likelihood, not maximum a posteriori) for each mode
        """
        out = dict()
        for mode in self._analyzer.get_mode_stats()['modes']:
            out[mode['index']] = {"parameters": mode['maximum']}
            self._set_model_parameters(mode['maximum'])
            out[mode['index']]["chi2"] = self._event.get_chi2()
            # The above line is not best performence, but easy solution.

        self._best_models_for_modes_MN = out

    def _parse_results(self):
        """
        Call the function that prints and saves results
        """
        if self._fit_method == "EMCEE":
            self._parse_results_EMCEE()
            if self._posterior_file_name is not None:
                self._save_posterior_EMCEE()
        elif self._fit_method == "MultiNest":
            self._parse_results_MultiNest()
        elif self._fit_method == "UltraNest":
            self._parse_results_UltraNest()
        else:
            raise ValueError('internal bug')

        # Below we close open files and remove temporary ones.
        if self._yaml_results:
            if self._yaml_results_file is not sys.stdout:
                self._yaml_results_file.close()
        if self._fit_method == "MultiNest":
            if self._MN_temporary_files:
                shutil.rmtree(self._kwargs_MultiNest['outputfiles_basename'],
                              ignore_errors=True)

    def _parse_results_EMCEE(self):
        """
        Print and save results from EMCEE fitting.

        This version works with EMCEE version 2.X and 3.0.
        """
        if self._yaml_results:
            lst = [mm.__version__, __version__]
            code_version = "MulensModel and script versions: {:}".format(lst)
            print(code_version, **self._yaml_kwargs)

        accept_rate = np.mean(self._sampler.acceptance_fraction)
        out = "Mean acceptance fraction: {0:.3f}".format(accept_rate)
        print(out)
        if self._yaml_results:
            print(out, **self._yaml_kwargs)

        auto_time = np.mean(self._sampler.get_autocorr_time(
            quiet=True, discard=self._fitting_parameters['n_burn']))
        out = "Mean autocorrelation time [steps]: {0:.1f}".format(auto_time)
        print(out)
        if self._yaml_results:
            print(out, **self._yaml_kwargs)
        self._extract_posterior_samples_EMCEE()

        if self._yaml_results and isinstance(self._fixed_parameters, dict):
            print("Fixed parameters:", **self._yaml_kwargs)
            for (key, value) in self._fixed_parameters.items():
                print("    {:} : {:}".format(key, value), **self._yaml_kwargs)

        print("Fitted parameters:")
        self._print_results(self._samples_flat)
        if self._yaml_results:
            print("Fitted parameters:", **self._yaml_kwargs)
            self._print_yaml_results(self._samples_flat)
        self._shift_t_0_in_samples()

        if self._return_fluxes:
            print("Fitted fluxes (source and blending):")
            blob_samples = self._get_fluxes_to_print_EMCEE()
            self._print_results(blob_samples, names='fluxes')
            if self._yaml_results:
                print("Fitted fluxes: # (source and blending)",
                      **self._yaml_kwargs)
                self._print_yaml_results(blob_samples, names='fluxes')

        self._print_best_model()
        if self._yaml_results:
            self._print_yaml_best_model()

        if self._shift_t_0 and self._yaml_results:
            print("Plots shift_t_0 : {:}".format(self._shift_t_0_val),
                  **self._yaml_kwargs)

    def _extract_posterior_samples_EMCEE(self):
        """
        set self._samples_flat and self._samples for EMCEE
        """
        n_burn = self._fitting_parameters['n_burn']
        self._samples = self._sampler.chain[:, n_burn:, :]
        for parameter in self._parse_posterior_abs:
            index = self._fit_parameters.index(parameter)
            self._samples[:, :, index] = np.fabs(self._samples[:, :, index])

        n_fit = self._n_fit_parameters
        self._samples_flat = self._samples.copy().reshape((-1, n_fit))
        if 'trace' not in self._plots:
            self._samples = None

    def _print_results(self, data, names="parameters", mode=None):
        """
        calculate and print median values and +- 1 sigma for given parameters
        """
        if names == "parameters":
            ids = self._fit_parameters
        elif names == "fluxes":
            if self._flux_names is None:
                self._flux_names = self._get_fluxes_names_to_print()
            ids = self._flux_names
        else:
            raise ValueError('internal bug')

        if self._fit_method == "EMCEE":
            results = self._get_weighted_percentile(data)
        elif self._fit_method in ["MultiNest", "UltraNest"]:
            if mode is None:
                weights = self._samples_flat_weights
            else:
                weights = self._samples_modes_flat_weights[mode]
            results = self._get_weighted_percentile(data, weights=weights)
        else:
            raise ValueError("internal bug")

        print(self._format_results(ids, results))

    def _format_results(self, ids, results, yaml=False, begin=""):
        """
        take a list of parameters and a list of results and
        return properly formatted string
        """
        text = ""
        for (parameter, results_) in zip(ids, results):
            format_ = "{:} : {:.5f} +{:.5f} -{:.5f}\n"
            if parameter != 'q':
                format_ = "{:} : {:.5f} +{:.5f} -{:.5f}\n"
                if yaml:
                    format_ = "{:} : [{:.5f}, +{:.5f}, -{:.5f}]\n"
            else:
                format_ = "{:} : {:.7f} +{:.7f} -{:.7f}\n"
                if yaml:
                    format_ = "{:} : [{:.7f}, +{:.7f}, -{:.7f}]\n"
            if parameter in self._parse_posterior_abs:
                parameter = "|{:}|".format(parameter)
            text += (begin + format_).format(parameter, *results_)
        return text[:-1]

    def _print_yaml_results(self, data, names="parameters", mode=None):
        """
        calculate and print in yaml format median values and +- 1 sigma
        for given parameters
        """
        begin = "    "
        if names == "parameters":
            ids = self._fit_parameters
        elif names == "fluxes":
            if self._flux_names is None:
                self._flux_names = self._get_fluxes_names_to_print()
            ids = self._flux_names
        else:
            raise ValueError('internal bug')

        if self._fit_method == "EMCEE":
            results = self._get_weighted_percentile(data)
        elif self._fit_method in ["MultiNest", "UltraNest"]:
            if mode is None:
                weights = self._samples_flat_weights
            else:
                weights = self._samples_modes_flat_weights[mode]
            results = self._get_weighted_percentile(data, weights=weights)
        else:
            raise ValueError("internal bug")

        print("# [median, sigma+, sigma-]", **self._yaml_kwargs)
        print(self._format_results(ids, results, yaml=True, begin=begin),
              **self._yaml_kwargs)

    def _get_fluxes_names_to_print(self):
        """
        get strings to be used as names of parameters to be printed
        """
        if self._n_fluxes_per_dataset == 2:
            s_or_b = ['s', 'b']
        elif self._n_fluxes_per_dataset == 3:
            s_or_b = ['s1', 's2', 'b']
        else:
            raise ValueError(
                'Internal error: ' + str(self._n_fluxes_per_dataset))
        n = self._n_fluxes_per_dataset
        flux_names = ['flux_{:}_{:}'.format(s_or_b[i % n], i // n+1)
                      for i in range(self._n_fluxes)]
        return flux_names

    def _get_weighted_percentile(
            self, data, fractions=[0.158655, 0.5, 0.841345], weights=None):
        """
        Calculate weighted percentile of the data. Data can be weighted or not.
        """
        if weights is None:
            kwargs = dict()
            if data.shape[0] > data.shape[1]:
                kwargs['axis'] = 0
            results = np.percentile(data, 100.*np.array(fractions), **kwargs)
        else:
            results = []
            for i in range(data.shape[1]):
                data_ = data[:, i]
                indexes = np.argsort(data_)
                weights_sorted = weights[indexes]
                weights_cumulative = np.cumsum(weights_sorted)
                weighted_quantiles = weights_cumulative - 0.5 * weights_sorted
                weighted_quantiles /= weights_cumulative[-1]
                results.append(
                    np.interp(fractions, weighted_quantiles, data_[indexes]))
            results = np.array(results).T

        out = []
        for i in range(results.shape[1]):
            median = results[1, i]
            out.append([median, results[2, i]-median, median-results[0, i]])

        return out

    def _shift_t_0_in_samples(self):
        """
        shift the values of t_0, t_0_1, and t_0_2:
        """
        if not self._shift_t_0:
            return

        for name in ['t_0', 't_0_1', 't_0_2']:
            if name in self._fit_parameters:
                index = self._fit_parameters.index(name)
                values = self._samples_flat[:, index]
                self._shift_t_0_val = int(np.mean(values))
                try:
                    self._samples_flat[:, index] -= self._shift_t_0_val
                    if 'trace' in self._plots:
                        self._samples[:, :, index] -= self._shift_t_0_val
                except TypeError:
                    fmt = ("Warning: extremely wide range of posterior {:}: "
                           "from {:} to {:}")
                    warnings.warn(
                        fmt.format(name, np.min(values), np.max(values)))
                    self._samples_flat[:, index] = values - self._shift_t_0_val
                    if 'trace' in self._plots:
                        self._samples[:, :, index] = (
                            self._samples[:, :, index] - self._shift_t_0_val)

            if self._fixed_parameters is not None:
                if name in self._fixed_parameters.keys():
                    self._shift_t_0_val = int(self._fixed_parameters[name])

    def _get_fluxes_to_print_EMCEE(self):
        """
        prepare values to be printed for EMCEE fitting
        """
        try:
            blobs = np.array(self._sampler.blobs)
        except Exception as exception:
            raise ValueError('There was some issue with blobs:\n' +
                             str(exception))
        blob_sampler = np.transpose(blobs, axes=(1, 0, 2))
        blob_samples = blob_sampler[:, self._fitting_parameters['n_burn']:, :]
        blob_samples = blob_samples.reshape((-1, self._n_fluxes))

        return blob_samples

    def _print_best_model(self):
        """
        print best model found
        """
        print("Best model:")
        if self._flat_priors:
            print("chi2: {:.4f}".format(-2. * self._best_model_ln_prob))
        else:
            self._ln_like(self._best_model_theta)
            print("chi2: {:.4f}".format(self._event.get_chi2()))
        print(*self._fit_parameters)
        print(*list(self._best_model_theta))
        if self._return_fluxes:
            print("Fluxes:")
            print(*list(self._best_model_fluxes))

    def _print_yaml_best_model(self, begin="", mode=None):
        """
        print in yaml format best model found
        """
        yaml_txt = begin + "Best model:\n"

        if mode is None:
            zip_1 = zip(self._fit_parameters, self._best_model_theta)
            zip_2 = zip(self._flux_names, self._best_model_fluxes)
            if self._flat_priors:
                chi2 = -2. * self._best_model_ln_prob
            else:
                self._ln_like(self._best_model_theta)
                chi2 = self._event.get_chi2()
        else:
            zip_1 = zip(self._fit_parameters,
                        mode['parameters'][:self._n_fit_parameters])
            zip_2 = zip(self._flux_names,
                        mode['parameters'][self._n_fit_parameters:])
            chi2 = mode['chi2']

        yaml_txt += (begin + "  chi2: {:.4f}\n").format(chi2)
        yaml_txt += begin + "  Parameters:\n"
        format_ = begin + "    {:}: {:}\n"
        for (parameter, results_) in zip_1:
            yaml_txt += format_.format(parameter, results_)

        if self._flux_names is None:
            self._flux_names = self._get_fluxes_names_to_print()
        yaml_txt += begin + "  Fluxes:\n"
        for (parameter, results_) in zip_2:
            yaml_txt += format_.format(parameter, results_)

        print(yaml_txt, end="", **self._yaml_kwargs)

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

    def _parse_results_MultiNest(self):
        """
        Parse results of MultiNest fitting
        """
        self._extract_posterior_samples_MultiNest()

        if self._kwargs_MultiNest['multimodal']:
            self._parse_results_MultiNest_multimodal()
        else:
            self._parse_results_MultiNest_singlemode()

        self._shift_t_0_in_samples()

        self._print_best_model()
        if self._yaml_results:
            self._print_yaml_best_model()

    def _parse_results_MultiNest_multimodal(self):
        """
        Print parameters and fluxes for each mode separately
        """
        out = "Number of modes found: {:}".format(self._n_modes)
        print(out)
        if self._yaml_results:
            print(out, **self._yaml_kwargs)
        if self._return_fluxes:
            print("Fitted parameters and fluxes (source and blending) "
                  "plus best model info:")
        else:
            print("Fitted parameters:")

        self._set_mode_probabilities()

        for i_mode in range(self._n_modes):
            self._parse_results_MultiNest_multimodal_one_mode(i_mode)

        print(" END OF MODES")
        if self._yaml_results:
            print("# END OF MODES", **self._yaml_kwargs)

    def _parse_results_MultiNest_multimodal_one_mode(self, i_mode):
        """
        Print parameters and fluxes for one of many modes
        """
        probability = self._mode_probabilities[i_mode]
        err = self._mode_probabilities_error[i_mode]
        accuracy = self._get_accuracy(err)
        samples = self._samples_modes_flat[i_mode]
        fluxes = self._samples_modes_flat_fluxes[i_mode]
        mode = self._best_models_for_modes_MN[i_mode]

        fmt = (" MODE {:} probability: {:." + str(accuracy) +
               "f} +- {:." + str(accuracy) + "f}")
        print(fmt.format(i_mode+1, probability, err))
        self._print_results(samples, mode=i_mode)

        if self._yaml_results:
            fmt = ("MODE {:}:\n  probability: {:}\n  "
                   "probability_sigma: {:}\n  Fitted parameters: ")
            print(fmt.format(i_mode+1, probability, err),
                  end="", **self._yaml_kwargs)
            self._print_yaml_results(samples, mode=i_mode)

        if self._return_fluxes:
            self._print_results(fluxes, names="fluxes")
            if self._yaml_results:
                print("  Fitted fluxes: # source and blending ", end="",
                      **self._yaml_kwargs)
                self._print_yaml_results(fluxes, names="fluxes")

        print("{:.4f}".format(mode['chi2']))
        print(*mode['parameters'][:self._n_fit_parameters])
        print(*mode['parameters'][self._n_fit_parameters:])
        if self._yaml_results:
            self._print_yaml_best_model(mode=mode, begin="  ")

    def _parse_results_MultiNest_singlemode(self):
        """
        Print results for MultiNest run in "single mode" setup.
        Note that this is different than
        _parse_results_MultiNest_multimodal_one_mode()
        """
        print("Fitted parameters:")
        self._print_results(self._samples_flat)

        if self._yaml_results:
            print("Fitted parameters:", **self._yaml_kwargs)
            self._print_yaml_results(self._samples_flat)

        if self._return_fluxes:
            print("Fitted fluxes (source and blending):")
            flux_samples = self._get_fluxes_to_print_MultiNest()
            self._print_results(flux_samples, names='fluxes')
            if self._yaml_results:
                print("Fitted fluxes: # (source and blending)",
                      **self._yaml_kwargs)
                self._print_yaml_results(flux_samples, names='fluxes')

    def _set_mode_probabilities(self):
        """
        Calculate probabilities of each mode and its uncertainty
        """
        modes = self._analyzer.get_mode_stats()['modes']
        key_lnZ = 'local log-evidence'
        modes_lnZ = np.array([mode[key_lnZ] for mode in modes])
        shift = (np.max(modes_lnZ) + np.min(modes_lnZ)) / 2.
        # We subtract shift for numerical stability.
        modes_lnZ -= shift
        modes_Z = np.exp(modes_lnZ)
        self._mode_probabilities = modes_Z / np.sum(modes_Z)
        relative_error = [mode[key_lnZ + ' error'] for mode in modes]
        # Error in np.sum(modes_Z) is ignored.
        self._mode_probabilities_error = (
            relative_error * self._mode_probabilities)

    def _get_accuracy(self, uncertainty):
        """
        Return int which says how to round given value to 2 significant digits
        """
        return 2-int(np.log10(uncertainty))

    def _extract_posterior_samples_MultiNest(self):
        """
        set self._samples_flat and self._samples_flat_weights for MultiNest
        """
        index = 2 + self._n_fit_parameters

        self._samples_flat = self._analyzer_data[:, 2:index]
        self._samples_flat_weights = self._analyzer_data[:, 0]

        if self._kwargs_MultiNest['multimodal']:
            self._samples_modes_flat = []
            self._samples_modes_flat_weights = []
            self._samples_modes_flat_fluxes = []
            n_begin = 0
            for n in self._MN_modes_indexes:
                n_end = n_begin + n
                weights = self._MN_samples_modes_all[n_begin:n_end, 0]
                samples = self._MN_samples_modes_all[n_begin:n_end, 2:index]
                fluxes_ = self._MN_samples_modes_all[n_begin:n_end, index:]
                self._samples_modes_flat.append(samples)
                self._samples_modes_flat_weights.append(weights)
                self._samples_modes_flat_fluxes.append(fluxes_)
                n_begin = n_end

    def _get_fluxes_to_print_MultiNest(self):
        """
        prepare flux values to be printed for MultiNest and Ultranest fitting
        """
        if self._fit_method == "MultiNest":
            index = 2 + self._n_fit_parameters
            data = self._analyzer_data
        elif self._fit_method == "UltraNest":
            index = self._n_fit_parameters
            data = self._result_UltraNest['weighted_samples']['points']

        return data[:, index:]

    def _parse_results_UltraNest(self):
        """
        Parse results of UltraNest fitting.
        Functions that print and save EMCEE results are also called here.
        """
        # re-weighted posterior samples:
        # self._samples_flat = self._result_UltraNest['samples'][:, :-2]
        # weighted samples from the posterior:
        weighted_samples = self._result_UltraNest['weighted_samples']
        index = self._n_fit_parameters
        self._samples_flat = weighted_samples['points'][:, :index]
        self._samples_flat_weights = weighted_samples['weights']
        self._sampler.print_results()

        max_like = self._result_UltraNest['maximum_likelihood']
        self._best_model_ln_prob = max_like['logl']
        self._best_model_theta = max_like['point'][:self._n_fit_parameters]
        self._best_model_fluxes = max_like['point'][self._n_fit_parameters:]
        self._parse_results_MultiNest_singlemode()

        self._shift_t_0_in_samples()
        self._print_best_model()
        if self._yaml_results:
            self._print_yaml_best_model()
            ln_ev = self._result_UltraNest['logz_single']
            ln_ev_err = self._result_UltraNest['logzerr_single']
            lns = "  ln_ev: [{:.5f}, +{:.5f}, -{:.5f}]"
            print(lns.format(ln_ev, ln_ev_err, ln_ev_err), **self._yaml_kwargs)

    def _write_residuals(self):
        """
        Write residuals to a file
        """
        if not self._residuals_output:
            return

        zip_ = zip(self._datasets, self._residuals_files, self._event.fits)
        for (dataset, name, fit) in zip_:
            if name == "-":
                continue
            (residuals, uncertainties) = fit.get_residuals(
                phot_fmt=dataset.input_fmt)
            data_out = np.array([dataset.time, residuals, uncertainties]).T
            np.savetxt(name, data_out)
            args = [1+self._datasets.index(dataset), name]
            print("Residuals for file number {:} saved in {:}".format(*args))

    def _make_plots(self):
        """
        make plots after fitting: best model, triangle plot, trace plot
        """
        if 'triangle' in self._plots:
            self._triangle_plot()
        if 'trace' in self._plots:
            self._trace_plot()
        if 'best model' in self._plots:
            self._best_model_plot()
            if 'interactive' in self._plots['best model']:
                self._make_interactive_plot()
        if 'trajectory' in self._plots:
            self._make_trajectory_plot()

    def _triangle_plot(self):
        """
        Make a triangle plot
        """
        self._reset_rcParams()

        n_bins = 40

        kwargs = {
            'bins': n_bins, 'show_titles': True, 'top_ticks': False,
            'quantiles': [0.15866, 0.5, 0.84134], 'verbose': False}

        data = self._get_samples_for_triangle_plot()
        labels = self._get_labels_for_triangle_plot()
        if data.shape[1] != len(labels):
            msg = "This should never happen: {:} {:}"
            raise ValueError(msg.format(data.shape, len(labels)))

        figure = corner.corner(data, labels=labels, **kwargs)

        self._save_figure(self._plots['triangle'].get('file'), figure=figure)

    def _reset_rcParams(self):
        """
        Reset matplotlib rcParams to their defaults
        """
        rcParams.update(rcParamsDefault)

    def _get_samples_for_triangle_plot(self):
        """
        Prepare samples that will be plotted on triangle plot
        """
        return self._samples_flat

    def _get_labels_for_triangle_plot(self):
        """
        provide list of labels to be used by triangle plot
        """
        return self._fit_parameters_latex

    def _save_figure(self, file_name, figure=None, dpi=None):
        """
        Save figure or display it
        """
        if file_name is None:
            plt.show()
        # XXX - does this work?
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
            if path.isfile(file_name):
                msg = "Existing file " + file_name + " will be overwritten"
                warnings.warn(msg)
            caller.savefig(file_name, **kwargs)
        plt.close()

    def _trace_plot(self):
        """
        Make a trace plot
        """
        self._reset_rcParams()

        alpha = 0.5
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.4
        figure_size = (7.5, 10.5)
        margins = {'left': 0.13, 'right': 0.97, 'top': 0.98, 'bottom': 0.05}

        n_grid = self._n_fit_parameters
        if len(self._fit_parameters_latex) != n_grid:
            msg = "This should never happen: {:} {:}"
            raise ValueError(
                msg.format(len(self._fit_parameters_latex), n_grid))
        grid = gridspec.GridSpec(n_grid, 1, hspace=0)

        plt.figure(figsize=figure_size)
        plt.subplots_adjust(**margins)
        x_vector = np.arange(self._samples.shape[1])

        for (i, latex_name) in enumerate(self._fit_parameters_latex):
            if i == 0:
                plt.subplot(grid[i])
                ax0 = plt.gca()
            else:
                plt.gcf().add_subplot(grid[i], sharex=ax0)
            plt.ylabel(latex_name)
            for j in range(self._samples.shape[0]):
                plt.plot(x_vector, self._samples[j, :, i], alpha=alpha)
            plt.xlim(0, self._samples.shape[1])
            plt.gca().tick_params(axis='both', which='both', direction='in',
                                  top=True, right=True)
            if i != self._n_fit_parameters - 1:
                plt.setp(plt.gca().get_xticklabels(), visible=False)
            plt.gca().set_prop_cycle(None)
        plt.xlabel('step count')

        self._save_figure(self._plots['trace'].get('file'))

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

        self._plot_models_for_best_model_plot(fluxes, kwargs_model)

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
        default_model = {'color': 'black', 'linewidth': 1.0, 'zorder': np.inf}

        kwargs_grid = {
            'nrows': 2, 'ncols': 1, 'height_ratios': [plot_size_ratio, 1],
            'hspace': hspace}

        (t_1, t_2) = self._get_time_limits_for_plot(tau, 'best model')
        (kwargs, xlim) = self._get_subtract_kwargs_and_xlim_for_best_model(t_1, t_2)

        kwargs_model = {
            't_start': t_1, 't_stop': t_2, **default_model, **kwargs}
        if self._model.n_sources != 1:
            kwargs_model['source_flux_ratio'] = self._datasets[0]
        if self._datasets[0].bandpass is not None:
            key = 'limb darkening u'
            if self._datasets[0].bandpass in self._model_parameters[key]:
                u = self._model_parameters[key][self._datasets[0].bandpass]
                kwargs_model['gamma'] = mm.Utils.u_to_gamma(u)

        kwargs_axes_1 = dict(
            axis='both', direction='in', bottom=True, top=True, left=True,
            right=True, labelbottom=False)
        kwargs_axes_2 = {**kwargs_axes_1, 'labelbottom': True}

        return (kwargs_grid, kwargs_model, kwargs, xlim, t_1, t_2,
                kwargs_axes_1, kwargs_axes_2)

    def _get_subtract_kwargs_and_xlim_for_best_model(self, t_1, t_2):
        """
        Should we subtract 2450000 or 2460000?
        """
        mean = (t_1 + t_2) / 2.
        if mean < 2460000.:
            kwargs = {'subtract_2450000': True, 'subtract_2460000': False}
            xlim = [t_1-2450000., t_2-2450000.]
        else:
            kwargs = {'subtract_2450000': False, 'subtract_2460000': True}
            xlim = [t_1-2460000., t_2-2460000.]

        return (kwargs, xlim)

    def _get_time_limits_for_plot(self, tau, plot_type):
        """
        find limits for the best model or trajectory plot
        """
        if 'time range' in self._plots[plot_type]:
            t_1 = self._plots[plot_type]['time range'][0]
            t_2 = self._plots[plot_type]['time range'][1]
            return (t_1, t_2)

        if self._model.n_sources == 1:
            t_1 = self._model.parameters.t_0
            t_2 = self._model.parameters.t_0
        elif self._model.n_sources == 2:
            if self._model.parameters.is_xallarap:
                t_1 = self._model.parameters.t_0
                t_2 = self._model.parameters.t_0
            else:
                t_1 = self._model.parameters.t_0_1
                t_2 = self._model.parameters.t_0_2
                if t_1 > t_2:
                    (t_1, t_2) = (t_2, t_1)
        else:
            raise ValueError('internal issue: ' + str(self._model.n_sources))

        t_1 -= tau * self._model.parameters.t_E
        t_2 += tau * self._model.parameters.t_E

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
            if np.sum(mask) == 0:
                continue
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

    def _plot_models_for_best_model_plot(self, fluxes, kwargs_model):
        """
        Plot best models: first ground-based (if needed, hence loop),
        then satellite ones (if needed).
        """
        for dataset in self._datasets:
            if dataset.ephemerides_file is None:
                self._model.plot_lc(source_flux=fluxes[0], blend_flux=fluxes[1], **kwargs_model)
                break

        for model in self._models_satellite:
            model.parameters.parameters = {**self._model.parameters.parameters}
            model.plot_lc(source_flux=fluxes[0], blend_flux=fluxes[1], **kwargs_model)

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
        (magnifications, labels, ylim, ax2) = self._second_Y_axis_settings()
        (A_range, ref_fluxes) = self._second_Y_axis_get_fluxes(ylim)
        out1, out2 = False, False
        if magnifications == "optimal":
            (magnifications, labels, out1) = self._second_Y_axis_optimal(
                ax2, *A_range)
            self._second_Y_axis_minor_ticks(ax2, magnifications, ref_fluxes)
        flux = ref_fluxes[0] * np.array(magnifications) + ref_fluxes[1]
        out2 = self._second_Y_axis_warnings(flux, labels, magnifications,
                                            *A_range)
        if out1 or out2:
            ax2.get_yaxis().set_visible(False)
            return

        ticks = mm.Utils.get_mag_from_flux(flux)
        try:  # matplotlib version 3.5 or later
            ax2.set_yticks(ticks=ticks, labels=labels)
        except Exception:  # matplotlib version 3.4.X or smaller
            ax2.set_yticks(ticks=ticks)
            ax2.set_yticklabels(labels=labels)

        ax2.set_ylim(ylim[0], ylim[1])

    def _second_Y_axis_settings(self):
        """
        Get and apply settings for the second Y axis
        """
        settings = self._plots['best model']["second Y scale"]
        magnifications = settings['magnifications']
        color = settings.get("color", "black")
        label = settings.get("label", "Magnification")
        labels = settings.get("labels")
        ylim = plt.ylim()

        ax2 = plt.gca().twinx()
        ax2.set_ylabel(label).set_color(color)
        ax2.spines['right'].set_color(color)
        ax2.tick_params(axis='y', direction="in", which="both", colors=color)

        return (magnifications, labels, ylim, ax2)

    def _second_Y_axis_get_fluxes(self, ylim):
        """
        Get fluxes and limiting magnification values for the second Y axis
        """
        flux_min = mm.Utils.get_flux_from_mag(ylim[0])
        flux_max = mm.Utils.get_flux_from_mag(ylim[1])
        (source_flux, blend_flux) = self._event.get_ref_fluxes()

        total_source_flux = sum(source_flux)
        A_min = (flux_min - blend_flux) / total_source_flux
        A_max = (flux_max - blend_flux) / total_source_flux

        return ([A_min, A_max], [total_source_flux, blend_flux])

    def _second_Y_axis_optimal(self, ax2, A_min, A_max):
        """
        Get optimal values of magnifications and labels
        """
        ax2.set_ylim(A_min, A_max)
        A_values = ax2.yaxis.get_ticklocs().round(7)
        A_values = A_values[(A_values >= max(1, A_min)) & (A_values < A_max)]
        if 1. not in A_values and A_min <= 1:
            A_values = np.insert(A_values, 0, 1.)
        is_integer = [mag.is_integer() for mag in A_values]
        if all(is_integer):
            labels = [f"{int(x):d}" for x in A_values]
            return (A_values, labels, False)

        fnum = np.array([str(x)[::-1].find(".") for x in A_values])
        labels = np.array([f"%0.{max(fnum)}f" % x for x in A_values])
        if max(fnum) >= 4 and len(fnum[fnum < 4]) < 3:
            msg = ("The computed magnifications for the second Y scale cover"
                   " a range too small to be shown: {:}")
            warnings.warn(msg.format(A_values))
            return (A_values, labels.tolist(), True)
        if max(fnum) >= 4:
            labels = np.array([f"{x:0.3f}" for x in A_values])

        return (A_values[fnum < 4], labels[fnum < 4].tolist(), False)

    def _second_Y_axis_minor_ticks(self, ax2, A_values, ref_fluxes):
        """
        Get minor ticks for magnification axis from matplotlib
        """
        ax2.minorticks_on()
        minor_ticks_A = np.array(ax2.yaxis.get_ticklocs(minor=True))
        minor_ticks_A = minor_ticks_A[~np.isin(minor_ticks_A, A_values)]

        minor_ticks_flux = ref_fluxes[0] * minor_ticks_A + ref_fluxes[1]
        minor_ticks_mag = mm.Utils.get_mag_from_flux(minor_ticks_flux)
        ax2.set_yticks(minor_ticks_mag, minor=True)

    def _second_Y_axis_warnings(self, flux, labels, A_values, A_min, A_max):
        """
        Issue warnings for negative flux or bad range of magnificaitons
        """
        if np.any(flux < 0.):
            mask = (flux > 0.)
            flux = flux[mask]
            labels = [l for (l, m) in zip(labels, mask) if m]
            msg = ("\n\n{:} label/s on the second Y scale will not be shown "
                   "because they correspond to negative flux which cannot "
                   "be translated to magnitudes.")
            warnings.warn(msg.format(np.sum(np.logical_not(mask))))

        if (np.min(A_values) < A_min or np.max(A_values) > A_max or
                np.any(flux < 0.)):
            msg = ("Provided magnifications for the second (i.e., right-hand "
                   "side) Y-axis scale are from {:} to {:},\nbut the range "
                   "of plotted magnifications is from {:} to {:}, hence, "
                   "the second scale is not plotted")
            args = [np.min(A_values), np.max(A_values), A_min, A_max]
            warnings.warn(msg.format(*args))
            return True

        return False

    def _make_trajectory_plot(self):
        """
        make and save plot of the best trajectory
        """
        dpi = 300
        tau = 1.5

        self._ln_like(self._best_model_theta)  # Sets all parameters to
        # the best model.

        self._reset_rcParams()

        if 'time range' in self._plots['trajectory']:
            t_range = self._get_time_limits_for_plot(tau, 'trajectory')
        else:
            t_range = self._get_time_limits_for_plot(tau, 'best model')
        kwargs = {'caustics': True, 't_range': t_range}

        self._model.plot_trajectory(**kwargs)

        plt.xlabel(r"$x [\theta_E]$")
        plt.ylabel(r"$y [\theta_E]$")
        plt.axis('equal')
        plt.gca().tick_params(axis='both', which='both', direction='in',
                              top=True, right=True)

        self._save_figure(self._plots['trajectory'].get('file'), dpi=dpi)

    def _make_interactive_plot(self):
        """
        plot best model and residuals interactively
        """
        scale = 0.5  # original size=(1920:1440)

        self._ln_like(self._best_model_theta)  # Sets all parameters to the best model.

        self._reset_rcParams()
        if 'rcParams' in self._plots['best model']:
            for (key, value) in self._plots['best model']['rcParams'].items():
                rcParams[key] = value

        kwargs_all = self._get_kwargs_for_best_model_plot()
        (ylim, ylim_residuals) = self._get_ylim_for_best_model_plot(*kwargs_all[4:6])
        (layout, kwargs_model, kwargs_interactive, kwargs) = \
            self._prepare_interactive_layout(scale, kwargs_all, ylim, ylim_residuals)

        (t_data_start, t_data_stop) = self._get_time_span_data()
        kwargs_model['t_start'] = t_data_start
        kwargs_model['t_stop'] = t_data_stop
        data_ref = self._event.data_ref
        (f_source_0, f_blend_0) = self._event.get_flux_for_dataset(data_ref)
        traces_lc = self._make_interactive_lc_traces(f_source_0, f_blend_0, **kwargs_model, **kwargs_interactive)
        self._interactive_fig = go.Figure(data=traces_lc, layout=layout)

        self._add_interactive_zero_trace(**kwargs_model, **kwargs_interactive)
        self._add_interactive_data_traces(kwargs_interactive, **kwargs)
        self._add_interactive_residuals_traces(kwargs_interactive, **kwargs_model)

        self._save_interactive_fig()

    def _prepare_interactive_layout(self, scale, kwargs_all, ylim, ylim_residuals):
        """Prepares the layout for the interactive plot."""
        kwargs_grid, kwargs_model, kwargs, xlim, t_1, t_2 = kwargs_all[:6]
        kwargs_axes_1, kwargs_axes_2 = kwargs_all[6:]
        kwargs_interactive = self._get_kwargs_for_plotly_plot(scale)

        layout = self._make_interactive_layout(
            ylim, ylim_residuals,
            **kwargs_grid,
            **kwargs_model,
            **kwargs_interactive
        )
        return layout, kwargs_model, kwargs_interactive, kwargs

    def _get_kwargs_for_plotly_plot(self, scale):
        """_
        setting kwargs for interactive plot
        """
        sizes = np.array([
            10.,  # markers data points
            4.,  # model line
            4.,  # residuals error thickens
            4.,  # residuals error width
            56.,  # font label
            4.,  # zero-line width in residuals
            4.,  # axes and legend border width
            15.,  # ticks len
            30.,  # font lagend
        ])
        sizes = sizes*scale
        colors = ['black', 'black', '#b9b9b9']  # This are: axes, font, and legend border.

        kwargs_interactive = dict(sizes=sizes, colors=colors, opacity=0.7, width=1920*scale,
                                  height=1440*scale, font='Old Standard TT, serif', paper_bgcolor='white')
        return kwargs_interactive

    def _make_interactive_layout(self, ylim, ylim_residuals, height_ratios, hspace, sizes, colors, opacity, width,
                                 height, font, paper_bgcolor, t_start, t_stop,
                                 subtract_2450000=None, subtract_2460000=None, **kwargs):
        """
        Creates plotly.graph_objects.Layout object analogues to best model plot
        """
        hsplit = height_ratios[1] / height_ratios[0]
        subtract = mm.utils.PlotUtils.find_subtract(subtract_2450000, subtract_2460000)

        xtitle = 'Time'
        if subtract > 0.:
            xtitle = xtitle+' - {:d}'.format(int(subtract))

        t_start = t_start - subtract
        t_stop = t_stop - subtract

        font_base = dict(family=font, size=sizes[4], color=colors[1])
        font_legend = dict(family=font, size=sizes[8])
        kwargs_ = dict(showgrid=False, ticks='inside', showline=True, ticklen=sizes[7],
                       tickwidth=sizes[6], linewidth=sizes[6], linecolor=colors[0], tickfont=font_base)
        kwargs_y = {'mirror': 'all', **kwargs_}
        kwargs_x = {'range': [t_start, t_stop], **kwargs_}
        layout = go.Layout(
            autosize=True, width=width, height=height, showlegend=True,
            legend=dict(
                x=1.02, y=.98, bgcolor=paper_bgcolor, bordercolor=colors[2], borderwidth=sizes[6], font=font_legend),
            paper_bgcolor=paper_bgcolor, plot_bgcolor=paper_bgcolor, font=font_base,
            yaxis=dict(title_text='Magnitude', domain=[hsplit+(hspace/2), 1], range=ylim, **kwargs_y),
            yaxis2=dict(title_text='Residuals', domain=[0, hsplit-(hspace/2)], anchor="x", range=ylim_residuals,
                        **kwargs_y),
            xaxis=dict(anchor="y", mirror='ticks', side='top', scaleanchor='x2', matches='x2', showticklabels=False,
                       **kwargs_x),
            xaxis2=dict(title_text=xtitle, anchor="y2", mirror='all', scaleanchor='x', matches='x', **kwargs_x)
            )
        return layout

    def _get_time_span_data(self):
        """
        Returning time span of datasets
        """
        t_min = np.zeros(len(self._datasets))
        t_max = np.zeros(len(self._datasets))
        for (i, data) in enumerate(self._datasets):
            t_min[i] = min(data.time)
            t_max[i] = max(data.time)

        return (min(t_min), max(t_max))

    def _make_interactive_lc_traces(self, f_source_0, f_blend_0, sizes, colors, opacity, width, height, font,
                                    paper_bgcolor, t_start, t_stop, name=None, dash='solid', subtract_2450000=None,
                                    subtract_2460000=None, gamma=None, bandpass=None, **kwargs):
        """
        Creates plotly.graph_objects.Scatter objects with model light curve
        """
        traces_lc = []
        subtract = mm.utils.PlotUtils.find_subtract(subtract_2450000, subtract_2460000)
        times = np.linspace(t_start, t_stop, num=5000) - subtract

        if isinstance(name, type(None)):
            showlegend = False
        else:
            showlegend = True

        for dataset in self._datasets:
            if dataset.ephemerides_file is None:
                lc = self._model.get_lc(
                    times=times, source_flux=f_source_0, blend_flux=f_blend_0, gamma=gamma, bandpass=bandpass)
                traces_lc.append(self._make_interactive_scatter_lc(
                    times, lc, name, showlegend, colors[1], sizes[1], dash))
                break

        traces_lc.extend(self._make_interactive_scatter_lc_satellite(
            traces_lc, times, f_source_0, f_blend_0, gamma, bandpass, colors, sizes, dash, subtract, showlegend))
        return traces_lc

    def _make_interactive_scatter_lc_satellite(
            self, traces, times, f_source_0, f_blend_0, gamma,
            bandpass, colors, sizes, dash, subtract, showlegend):
        """Generates Plotly Scatter traces for the light-curve satellite models."""

        for (i, model) in enumerate(self._models_satellite):
            name = self._event.datasets[i].plot_properties['label']
            model.parameters.parameters = {**self._model.parameters.parameters}
            lc = self._model.get_lc(times=times, source_flux=f_source_0, blend_flux=f_blend_0,
                                    gamma=gamma, bandpass=bandpass)
            times = times - subtract
            trace = self._make_interactive_scatter_lc(
                times, lc, name, showlegend, colors[1], sizes[1], dash)
            traces.append(trace)
        return traces

    def _make_interactive_scatter_lc(
            self, times, lc, name,
            showlegend, color, size, dash):
        """Creates a Plotly Scatter trace for the light curve."""

        return go.Scatter(x=times, y=lc, name=name, showlegend=showlegend, mode='lines',
                          line=dict(color=color, width=size, dash=dash),
                          xaxis="x", yaxis="y")

    def _add_interactive_zero_trace(self, t_start, t_stop, colors, sizes,
                                    subtract_2450000=False, subtract_2460000=False, **kwargs):
        """
        Creates plotly.graph_objects.Scatter object for line y=0 in
        residuals plot
        """
        subtract = mm.utils.PlotUtils.find_subtract(subtract_2450000, subtract_2460000)
        times = np.linspace(t_start, t_stop, num=2000)
        line = np.zeros(len(times))
        trace_0 = go.Scatter(x=times-subtract, y=line, mode='lines',
                             line=dict(color=colors[0], width=sizes[5], dash='dash'),
                             xaxis="x2", yaxis="y2", showlegend=False,)
        self._interactive_fig.add_trace(trace_0)

    def _add_interactive_data_traces(self, kwargs_interactive, phot_fmt='mag', data_ref=None, show_errorbars=True,
                                     show_bad=None, subtract_2450000=False, subtract_2460000=False, **kwargs):
        """
        Creates plotly.graph_objects.Scatter object for observation points
        per each data set.
        """
        self._event._set_default_colors()
        if self._event.fits is None:
            self._event.get_chi2()

        if data_ref is None:
            data_ref = self._event.data_ref

        subtract = mm.utils.PlotUtils.find_subtract(subtract_2450000, subtract_2460000)
        (f_source_0, f_blend_0) = self._event.get_flux_for_dataset(data_ref)

        traces_data = []
        for (dataset_index, data) in enumerate(self._datasets):
            # Scale the data flux
            (flux, err_flux) = self._event.fits[dataset_index].scale_fluxes(f_source_0, f_blend_0)
            (y_value, y_err) = mm.utils.PlotUtils.get_y_value_y_err(phot_fmt, flux, err_flux)
            times = data.time - subtract
            trace_data = self._make_one_interactive_data_trace(
                dataset_index, times, y_value, y_err, xaxis='x', yaxis='y', showlegend=True,
                show_errorbars=show_errorbars, show_bad=show_bad, **kwargs_interactive)
            traces_data.extend(trace_data)

        for trace in traces_data:
            self._interactive_fig.add_trace(trace)

    def _make_one_interactive_data_trace(self, dataset_index, times, y_value, y_err, xaxis, yaxis, showlegend,
                                         colors, sizes, opacity, show_errorbars=None, show_bad=None, **kwargs):
        """
        Creates plotly.graph_objects.Scatter object with data points form a given data set.
        """
        trace_data = []
        dataset, show_errorbars, show_bad = self._get_interactive_dataset(dataset_index)

        trace_data_good = self._make_interactive_good_data_trace(
            dataset, times, y_value, y_err, opacity, sizes, xaxis, yaxis, showlegend, show_errorbars)
        trace_data.append(trace_data_good)

        if show_bad:
            trace_data_bad = self._make_interactive_bad_data_trace(
                dataset, times, y_value, y_err, opacity, sizes, xaxis, yaxis, showlegend)
            trace_data.append(trace_data_bad)

        return trace_data

    def _get_interactive_dataset(self, dataset_index):
        """Get dataset properties for interactive plot settings."""
        dataset = self._event.datasets[dataset_index]
        show_errorbars = dataset.plot_properties.get('show_errorbars', True)
        show_bad = dataset.plot_properties.get('show_bad', False)
        return (dataset, show_errorbars, show_bad)

    def _make_interactive_good_data_trace(
            self, dataset, times, y_value, y_err, opacity,
            sizes, xaxis, yaxis, showlegend, show_errorbars):
        """Creates a single plotly.graph_objects.Scatter object for the good data points."""
        times_good = times[dataset.good]
        y_good = y_value[dataset.good]
        y_err_good = y_err[dataset.good]
        return self._make_interactive_data_trace(
            times_good, y_good, y_err_good, dataset,
            opacity, sizes, xaxis, yaxis, showlegend, show_errorbars)

    def _make_interactive_data_trace(self, x, y, y_err, dataset, opacity, sizes, xaxis, yaxis,
                                     showlegend, show_errorbars, color_override=None, error_visible=True):
        """Creates single plotly.graph_objects.Scatter object for good or bad data."""
        color = color_override if color_override else dataset.plot_properties['color']
        error_y = dict(type='data', array=y_err, visible=error_visible, thickness=sizes[2], width=sizes[3])
        marker = dict(color=color, size=sizes[0], line=dict(color=color, width=1))
        return go.Scatter(x=x, y=y, opacity=opacity, name=dataset.plot_properties['label'], mode='markers',
                          showlegend=showlegend, error_y=error_y, marker=marker, xaxis=xaxis, yaxis=yaxis)

    def _make_interactive_bad_data_trace(self, dataset, times, y_value, y_err, opacity, sizes,
                                         xaxis, yaxis, showlegend):
        """Creates a single plotly.graph_objects.Scatter object for the bad data points."""
        times_bad = times[dataset.bad]
        y_bad = y_value[dataset.bad]
        y_err_bad = y_err[dataset.bad]
        return self._make_interactive_data_trace(
            times_bad, y_bad, y_err_bad, dataset,
            opacity, sizes, xaxis, yaxis, showlegend,
            show_errorbars=False,
            color_override='black',
            error_visible=False
        )

    def _add_interactive_residuals_traces(self,  kwargs_interactive, phot_fmt='mag', data_ref=None, show_errorbars=True,
                                          show_bad=None, subtract_2450000=False, subtract_2460000=False, **kwargs,):
        """
        Creates plotly.graph_objects.Scatter object for residuals points
        per each data set.
        """
        traces_residuals = []
        self._event._set_default_colors()  # For each dataset
        if self._event.fits is None:
            self._event.get_chi2()

        if data_ref is None:
            data_ref = self._event.data_ref

        subtract = mm.utils.PlotUtils.find_subtract(subtract_2450000, subtract_2460000)

        # Get fluxes for the reference dataset
        (f_source_0, f_blend_0) = self._event.get_flux_for_dataset(data_ref)
        kwargs_residuals = {'phot_fmt': 'scaled', 'bad': False,
                            'source_flux': f_source_0, 'blend_flux': f_blend_0}
        if show_bad:
            kwargs_residuals['bad'] = True

        for (dataset_index, data) in enumerate(self._datasets):
            (y_value, y_err) = self._event.fits[dataset_index].get_residuals(**kwargs_residuals)
            times = data.time-subtract
            trace_residuals = self._make_one_interactive_data_trace(
                dataset_index, times, y_value, y_err, xaxis='x2', yaxis='y2',
                showlegend=False, show_errorbars=show_errorbars, show_bad=show_bad, **kwargs_interactive)
            traces_residuals.extend(trace_residuals)

        for trace in traces_residuals:
            self._interactive_fig.add_trace(trace)

    def _save_interactive_fig(self):
        """
        Saving interactive figure
        """
        file_ = self._plots['best model']['interactive']
        if path.exists(file_):
            if path.isfile(file_):
                msg = "Existing file " + file_ + " will be overwritten"
                warnings.warn(msg)
        self._interactive_fig.write_html(file_, full_html=True)


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
