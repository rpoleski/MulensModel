from astropy import units as u
import numpy as np
import warnings

from MulensModel.uniformcausticsampling import UniformCausticSampling


# For definition of class ModelParameters see below.

# Different parameter sets. Any parameters that may be given as
# 'basic' should be a list. Parameters that may be 'optional' should
# be a list of length 2. The second item will only be printed if the
# effect is included in the 'optional' list (see _get_effect_strings()).

_valid_parameters = {
    'point lens': ['t_0, u_0, t_E'],
    'point lens alt': 'alternate: t_eff may be substituted for u_0 or t_E',
    'binary lens': ['s, q, alpha'],
    'binary lens alt':
        'alternate: ' +
        '(x_caustic_in, x_caustic_out, t_caustic_in, t_caustic_out) ' +
        'may be substituted for (t_0, u_0, t_E, alpha)',
    'external_mass_sheet': ['convergence_K', 'shear_G', 'alpha'],
    'binary_lens_shear': ['convergence_K', 'shear_G'],
    'finite source': ['rho', '(for finite source effects)'],
    'finite source alt': 'alternate: t_star may be substituted for t_E or rho',
    'parallax': ['(pi_E_N, pi_E_E) OR pi_E', '(for parallax)'],
    'parallax opt': [
        't_0_par',
        'may also be specified for parallax models. Defaults to t_0.'],
    'lens orbital motion': ['dalpha_dt, ds_dt', '(for orbital motion)'],
    'lens orbital motion opt': [
        't_0_kep',
        'may also be specified for orbital motion models. Defaults to t_0.']}


def _get_effect_strings(*args):
    """
    Given *args[0], figure out which parameters should be printed.

    'basic' = fundamental parameters of the model or effect described
        by args[0]
    'additional' = any additional fundamental parameters (an extension
        of 'basic')
    'alternate' = possible substitutions for the fundamental parameters
    'optional' = parameters that may also be specified for the given
        model type.

    e.g. 'FSPL' returns basic = [t_0, u_0, tE], additional = [rho, s,
    q, alpha], alternate = [t_eff, t_star], and optional = [pi_E or
    pi_E_N, pi_E_E]
    """
    basic = None
    additional = []
    alternate = []
    optional = []

    args_0 = args[0].lower().replace(" ", "")

    # number of lenses
    if args_0 == 'pointlens' or args_0[2:4] == 'pl':
        basic = 'point lens'
        alternate.append('point lens alt')

    if args_0 == 'binarylens':
        basic = 'binary lens'

    if args_0[2:4] == 'bl':
        basic = 'point lens'
        additional.append('binary lens')
        alternate.append('point lens alt')
        alternate.append('binary lens alt')

    # Effects
    if args_0 == 'finitesource':
        basic = 'finite source'
        alternate.append('finite source alt')

    if args_0[0:2] == 'fs':
        additional.append('finite source')
        alternate.append('finite source alt')

    if args_0 == 'parallax':
        basic = 'parallax'
        optional.append('parallax opt')

    if len(args[0]) == 4:
        optional.append('parallax')
        optional.append('parallax opt')

    if args[0].lower() == 'lens orbital motion':
        basic = 'lens orbital motion'
        optional.append('lens orbital motion opt')

    if len(args[0]) == 4 and args[0][2:4].lower() == 'bl':
        optional.append('lens orbital motion')
        optional.append('lens orbital motion opt')

    return {
        'basic': basic, 'additional': additional, 'alternate': alternate,
        'optional': optional}


def _print_parameters(header, components):
    """
    Prints the given parameter information under the requested header.

    Arguments:
        header: *str*

        components: *dictionary*
            This should be created using _get_effect_strings()
    """
    print(header)
    if components['basic'] is not None:
        parameters_list = 'basic: {0}'.format(
            _valid_parameters[components['basic']][0])
        if len(components['additional']) > 0:
            for item in components['additional']:
                parameters_list += ', {0}'.format(_valid_parameters[item][0])
        print('{0}'.format(parameters_list))

    if len(components['alternate']) > 0:
        for item in components['alternate']:
            print('{0}'.format(_valid_parameters[item]))

    if len(components['optional']) > 0:
        for item in components['optional']:
            print('optional: {0} {1}'.format(
                _valid_parameters[item][0], _valid_parameters[item][1]))


def _print_all():
    """
    Give the user general information about common models and effects.
    """
    print('------------------------')
    print('Some common model types:')
    print('------------------------')
    _print_parameters('PSPL: ', _get_effect_strings('PSPL'))
    _print_parameters('----\nFSBL: ', _get_effect_strings('FSBL'))
    print('-----------------')
    print('Optional Effects:')
    print('-----------------')
    _print_parameters('finite source: ', _get_effect_strings('finite source'))
    _print_parameters('---------\nparallax: ', _get_effect_strings('parallax'))
    _print_parameters(
        '---------\nlens orbital motion: ',
        _get_effect_strings('lens orbital motion'))
    print('-----------------')
    print('All Options: (call using which_parameters([option]) )')
    print('-----------------')
    print("Model types: 'PSPL', 'FSPL', 'PSBL', 'FSBL'")
    print("Effects: 'point lens', 'binary lens', 'finite source', " +
          "'parallax', 'lens orbital motion'")


def which_parameters(*args):
    """
    Prints information on valid parameter combinations that can be
    used to define a model or a particular effect. May be called with
    no arguments (returns information on many types of models) or with
    one argument referring to a specific model (e.g. PSPL) or effect
    (e.g. parallax).

    Valid arguments: *str*
        Model types: 'PSPL', 'FSPL', 'PSBL', 'FSBL'

        Effects: 'point lens', 'binary lens', 'finite source',
        'parallax', 'lens orbital motion'
    """
    warnings.warn(
        "Warning: function which_parameters() does not show binary source " +
        "parameters!",
        RuntimeWarning)
    if len(args) == 0:
        _print_all()
    else:
        components = _get_effect_strings(*args)
        header = '---------\n{0} parameters:'.format(args[0])
        _print_parameters(header, components)


class ModelParameters(object):
    """
    A class for the basic microlensing model parameters (t_0, u_0,
    t_E, s, q, alpha, etc.). Can handle point lens or binary lens.
    The pi_E assumes NE coordinates (Parallel, Perpendicular
    coordinates are not supported).

    Arguments :
        parameters: *dictionary*
            A dictionary of parameters and their values. See
            :py:func:`which_parameters()` for valid parameter combinations.

    Attributes :
        parameters: *dictionary*
            A dictionary of parameters and their values. Do not use it to
            change parameter values, instead use e.g.:
            ``model_parameters.u_0 = 0.1`` or
            ``setattr(model_parameters, 'u_0', 0.1)``.

    Example:
        Define a point lens model:
            ``params = ModelParameters({'t_0': 2450000., 'u_0': 0.3,
            't_E': 35.})``

        Then you can print the parameters:
            ``print(params)``

    """

    def __init__(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError(
                'ModelParameters must be initialized with dict ' +
                "as a parameter\ne.g., ModelParameters({'t_0': " +
                "2456789.0, 'u_0': 0.123, 't_E': 23.45})")

        self._count_sources(parameters.keys())
        self._count_lenses(parameters.keys())
        self._set_type(parameters.keys())
        self._check_types('alpha' in parameters.keys())

        if self.n_sources == 1:
            self._check_valid_combination_1_source(parameters.keys())
            if self._type['Cassan08']:
                self._uniform_caustic = None
                self._standard_parameters = None
        elif self.n_sources == 2:
            self._check_valid_combination_2_sources(parameters.keys())
            if 't_E' not in parameters.keys():
                raise KeyError('Currently, the binary source calculations ' +
                               'require t_E to be directly defined, i.e., ' +
                               'has to be the same for both sources.')
            (params_1, params_2) = self._divide_parameters(parameters)
            try:
                self._source_1_parameters = ModelParameters(params_1)
            except Exception:
                print("ERROR IN ITIALIZING SOURCE 1")
                raise
            try:
                self._source_2_parameters = ModelParameters(params_2)
            except Exception:
                print("ERROR IN ITIALIZING SOURCE 2")
                raise
            # The block above forces checks from "== 1" block above to be
            # run on each source parameters separately.
        else:
            raise ValueError('wrong number of sources')
        self._set_parameters(parameters)

    def _count_sources(self, keys):
        """How many sources there are?"""
        binary_params = ['t_0_1', 't_0_2', 'u_0_1', 'u_0_2', 'rho_1', 'rho_2',
                         't_star_1', 't_star_2']
        common = set(binary_params).intersection(set(keys))
        if len(common) == 0:
            self._n_sources = 1
        elif len(common) == 1:
            raise ValueError('Wrong parameters - the only binary source ' +
                             'parameter is {:}'.format(common))
        else:
            common_no_1_2 = {param[:-2] for param in common}
            condition_1 = (len(common_no_1_2) == len(common))
            condition_2 = not (
                'rho' in common_no_1_2 and 't_star' in common_no_1_2)
            if condition_1 and condition_2:
                raise ValueError(
                    'Given binary source parameters do not allow defining ' +
                    'the Model: {:}'.format(common))
            self._n_sources = 2

    def _count_lenses(self, keys):
        """How many lenses there are?"""
        self._n_lenses = 1
        if 's' in keys or 'q' in keys:
            self._n_lenses = 2
        # Both standard and Cassen08 parameterizations require s and q

    def _set_type(self, keys):
        """
        sets self._type property, which indicates what type of a model we have
        """
        types = ['finite source', 'parallax', 'Cassan08',
                 'lens 2-parameter orbital motion', 'mass sheet', 'xallarap']
        out = {type_: False for type_ in types}

        temp = {
            'finite source': 'rho t_star rho_1 rho_2 t_star_1 t_star_2',
            'parallax': 'pi_E_N pi_E_E pi_E',
            'Cassan08':
                'x_caustic_in x_caustic_out t_caustic_in t_caustic_out',
            'lens 2-parameter orbital motion': 'dalpha_dt ds_dt',
            'mass sheet': 'convergence_K shear_G',
            'xallarap': ('xi_period xi_semimajor_axis xi_inclination '
                         'xi_Omega_node xi_argument_of_latitude_reference')}

        parameter_to_type = dict()
        for (key, values) in temp.items():
            for value in values.split():
                parameter_to_type[value] = key

        for key in keys:
            if key in parameter_to_type:
                out[parameter_to_type[key]] = True

        self._type = out

    def _check_types(self, alpha_defined):
        """
        Check if self._type values make sense
        """
        n_lenses = self._n_lenses
        n_sources = self._n_sources

        # Lens orbital motion requires binary lens:
        if self._type['lens 2-parameter orbital motion'] and n_lenses == 1:
            raise KeyError('Orbital motion of the lens requires two lens '
                           'components but only one was provided.')

        self._check_types_for_Cassan08()

        if alpha_defined:
            if self._n_lenses == 1 and not self._type['mass sheet']:
                raise KeyError(
                    'You defined alpha for single lens model '
                    'without external mass sheet. This is not allowed.')

        if n_sources > 1 and self._type['xallarap']:
            raise NotImplementedError('We have not yet implemented xallarap '
                                      'and multiple luminous sources')

    def _check_valid_combination_1_source(self, keys):
        """
        Check that the user hasn't over-defined the ModelParameters.
        """
        # Make sure that there are no unwanted keys
        allowed_keys = set((
            't_0 u_0 t_E t_eff rho t_star pi_E pi_E_N pi_E_E t_0_par '
            's q alpha dalpha_dt ds_dt t_0_kep convergence_K shear_G '
            't_0_1 t_0_2 u_0_1 u_0_2 rho_1 rho_2 t_star_1 t_star_2 '
            'x_caustic_in x_caustic_out t_caustic_in t_caustic_out '
            'xi_period xi_semimajor_axis xi_inclination xi_Omega_node '
            'xi_argument_of_latitude_reference xi_eccentricity '
            'xi_omega_periapsis t_0_xi').split())
        difference = set(keys) - allowed_keys
        if len(difference) > 0:
            derived_1 = ['gamma', 'gamma_perp', 'gamma_parallel']
            if set(keys).intersection(derived_1):
                msg = ('You cannot set gamma, gamma_perp, ' +
                       'or gamma_parallel. These are derived parameters. ' +
                       'You can set ds_dt and dalpha_dt instead.\n')
            else:
                msg = ""
            msg += 'Unrecognized parameters: {:}'.format(difference)
            raise KeyError(msg)

        if self._type['Cassan08']:
            self._check_valid_combination_1_source_Cassan08(keys)
        else:
            self._check_valid_combination_1_source_standard(keys)

    def _check_types_for_Cassan08(self):
        """
        Check if Cassan08 is used and if so, then make sure that
        the trajectory is rectilinear and there is only one source.
        """
        if not self._type['Cassan08']:
            return

        types = ['parallax', 'xallarap', 'lens 2-parameter orbital motion']
        for type_ in types:
            if self._type[type_]:
                raise NotImplementedError(
                    'Currently we do not allow Cassan (2008) '
                    'parameterization of binary lens and ' + type_)

        if self._n_sources > 1:
            raise NotImplementedError(
                "Cassan (2008) parameterization doesn't work for "
                "multi sources models")

    def _divide_parameters(self, parameters):
        """
        Divide an input dict into 2 - each source separately.
        Some of the parameters are copied to both dicts.
        """
        separate_parameters = (
            't_0_1 t_0_2 u_0_1 u_0_2 rho_1 rho_2 t_star_1 t_star_2'.split())
        parameters_1 = {}
        parameters_2 = {}
        for (key, value) in parameters.items():
            if key in separate_parameters:
                if key[-2:] == "_1":
                    parameters_1[key[:-2]] = value
                elif key[-2:] == "_2":
                    parameters_2[key[:-2]] = value
                else:
                    raise ValueError('unexpected error')
            else:
                parameters_1[key] = value
                parameters_2[key] = value
        return (parameters_1, parameters_2)

    def __repr__(self):
        """A nice way to represent a ModelParameters object as a string"""
        keys = self._get_keys_for_repr()
        formats = self._get_formats_dict_for_repr()
        ordered_keys = self._get_ordered_keys_for_repr()

        variables = ''
        values = ''
        for key in ordered_keys:
            if key not in keys:
                continue
            (full_name, value) = self._get_values_for_repr(formats[key], key)
            (fmt_1, fmt_2) = self._get_formats_for_repr(formats[key],
                                                        full_name)
            variables += fmt_1.format(full_name)
            values += fmt_2.format(value)

        return '{0}\n{1}'.format(variables, values)

    def _get_keys_for_repr(self):
        """
        get all the keys that will be printed
        """
        keys = set(self.parameters.keys())
        if 'pi_E' in keys:
            keys.remove('pi_E')
            keys |= {'pi_E_E', 'pi_E_N'}

        if 'pi_E_E' in keys or 'pi_E_N' in keys:
            keys |= {'t_0_par'}

        if 'ds_dt' in keys or 'dalpha_dt' in keys:
            keys |= {'t_0_kep'}

        if self.is_xallarap:
            keys |= {'t_0_xi'}

        return keys

    def _get_formats_dict_for_repr(self):
        """
        define formats that define how to print the numbers
        """
        # Below we define dict of dicts. Key of inner ones: 'width',
        # 'precision', and optional: 'unit' and 'name'.
        formats = {
            't_0': {'width': 13, 'precision': 5, 'unit': 'HJD'},
            'u_0': {'width': 9, 'precision': 6},
            't_eff': {'width': 10, 'precision': 6, 'unit': 'd'},
            't_E': {'width': 10, 'precision': 4, 'unit': 'd'},
            'rho': {'width': 7, 'precision': 5},
            't_star': {'width': 13, 'precision': 6, 'unit': 'd'},
            'pi_E_N': {'width': 9, 'precision': 5},
            'pi_E_E': {'width': 9, 'precision': 5},
            't_0_par': {'width': 13, 'precision': 5, 'unit': 'HJD'},
            's': {'width': 9, 'precision': 5},
            'q': {'width': 12, 'precision': 8},
            'alpha': {'width': 11, 'precision': 5, 'unit': 'deg'},
            'convergence_K': {'width': 12, 'precision': 8},
            'shear_G': {'width': 12, 'precision': 8},
            'ds_dt': {
                'width': 11, 'precision': 5, 'unit': '/yr', 'name': 'ds/dt'},
            'dalpha_dt': {
                'width': 18, 'precision': 5, 'unit': 'deg/yr',
                'name': 'dalpha/dt'},
            't_0_kep': {'width': 13, 'precision': 5, 'unit': 'HJD'},
            'x_caustic_in': {'width': 13, 'precision': 7},
            'x_caustic_out': {'width': 13, 'precision': 7},
            't_caustic_in': {'width': 13, 'precision': 5, 'unit': 'HJD'},
            't_caustic_out': {'width': 13, 'precision': 5, 'unit': 'HJD'},
            'xi_period': {'width': 10, 'precision': 4,
                          'unit': 'd', 'name': 'xallarap period'},
            'xi_semimajor_axis': {'width': 9, 'precision': 6,
                                  'name': 'xallarap semimajor axis'},
            'xi_inclination': {'width': 11, 'precision': 5, 'unit': 'deg',
                               'name': 'xallarap inclination'},
            'xi_Omega_node': {'width': 11, 'precision': 5, 'unit': 'deg',
                              'name': 'xallarap Omega node'},
            'xi_argument_of_latitude_reference': {
                'width': 11, 'precision': 5, 'unit': 'deg',
                'name': 'xallarap argument of latitude reference'},
            'xi_eccentricity': {'width': 8, 'precision': 6,
                                'name': 'xallarap eccentricity'},
            'xi_omega_periapsis': {'width': 11, 'precision': 5, 'unit': 'deg',
                                   'name': 'xallarap omega periapsis'},
            't_0_xi': {'width': 13, 'precision': 5, 'unit': 'HJD'},
        }
        # Add binary source parameters with the same settings.
        binary_source_keys = ['t_0_1', 't_0_2', 'u_0_1', 'u_0_2',
                              'rho_1', 'rho_2', 't_star_1', 't_star_2']
        for key in binary_source_keys:
            form = formats[key[:-2]]
            formats[key] = {'width': form['width'],
                            'precision': form['precision']}
            if 'unit' in form:
                formats[key]['unit'] = form['unit']
            if 'name' in form:
                raise KeyError('internal issue: {:}'.format(key))

        return formats

    def _get_ordered_keys_for_repr(self):
        """
        define the default order of parameters
        """
        ordered_keys = [
            't_0', 't_0_1', 't_0_2', 'u_0', 'u_0_1', 'u_0_2', 't_eff', 't_E',
            'rho', 'rho_1', 'rho_2', 't_star', 't_star_1', 't_star_2',
            'pi_E_N', 'pi_E_E', 't_0_par', 's', 'q', 'alpha',
            'convergence_K', 'shear_G', 'ds_dt', 'dalpha_dt', 't_0_kep',
            'x_caustic_in', 'x_caustic_out', 't_caustic_in', 't_caustic_out',
            'xi_period', 'xi_semimajor_axis', 'xi_inclination',
            'xi_Omega_node', 'xi_argument_of_latitude_reference',
            'xi_eccentricity', 'xi_omega_periapsis', 't_0_xi'
        ]
        return ordered_keys

    def _get_values_for_repr(self, form, key):
        """
        Get full name of the parameter and its value (float)
        to be used by __rerp__().
        """
        full_name = form.get('name', key)
        if 'unit' in form:
            full_name += " ({:})".format(form['unit'])

        value = getattr(self, key)
        if isinstance(value, u.Quantity):
            value = value.value

        return (full_name, value)

    def _get_formats_for_repr(self, form, full_name):
        """
        Extract formats to be used by __repr__().
        """
        fmt_1 = '{:>' + str(max([form['width'], len(full_name)]))
        fmt_2 = fmt_1 + '.' + str(form['precision']) + 'f} '
        fmt_1 += '} '
        return (fmt_1, fmt_2)

    def _check_valid_combination_2_sources(self, keys):
        """
        make sure that there is no conflict between t_0 and t_0_1 etc.
        """
        binary_params = (
            't_0_1 t_0_2 u_0_1 u_0_2 rho_1 rho_2 t_star_1 t_star_2'.split())
        for parameter in binary_params:
            if (parameter in keys) and (parameter[:-2] in keys):
                raise ValueError('You cannot set {:} and {:}'.format(
                                 parameter, parameter[:-2]))

    def _check_valid_combination_1_source_standard(self, keys):
        """
        Here we check parameters for non-Cassan08 parameterization.
        """
        self._check_valid_combination_1_source_t_0_u_0(keys)
        self._check_valid_combination_1_source_t_E(keys)
        self._check_valid_combination_1_source_parallax(keys)
        self._check_valid_combination_1_source_mass_sheet(keys)
        self._check_valid_combination_1_source_binary_lens(keys)
        self._check_valid_combination_1_source_xallarap(keys)

    def _check_valid_combination_1_source_t_0_u_0(self, keys):
        """
        Make sure that t_0 and u_0 are defined.
        """
        if 't_0' not in keys:
            raise KeyError('t_0 must be defined')

        if ('u_0' not in keys) and ('t_eff' not in keys):
            raise KeyError('not enough information to calculate u_0')

    def _check_valid_combination_1_source_t_E(self, keys):
        """
        Make sure that t_E is defined and that it's not overdefined.
        """
        if (('t_E' not in keys) and
                (('u_0' not in keys) or ('t_eff' not in keys)) and
                (('rho' not in keys) or ('t_star' not in keys))):
            raise KeyError('not enough information to calculate t_E')

        if (('rho' in keys) and ('t_star' in keys) and ('u_0' in keys) and
                ('t_eff' in keys)):
            raise KeyError('You cannot define rho, t_star, u_0, and t_eff')

        if ('t_E' in keys) and ('rho' in keys) and ('t_star' in keys):
            raise KeyError('Only 1 or 2 of (t_E, rho, t_star) may be defined.')

        if ('t_E' in keys) and ('u_0' in keys) and ('t_eff' in keys):
            raise KeyError('Only 1 or 2 of (u_0, t_E, t_eff) may be defined.')

    def _check_valid_combination_1_source_parallax(self, keys):
        """
        Here we check parallax parameters for non-Cassan08 parameterization.
        """
        # Parallax is either pi_E or (pi_E_N, pi_E_E)
        if 'pi_E' in keys and ('pi_E_N' in keys or 'pi_E_E' in keys):
            raise KeyError(
                'Parallax may be defined EITHER by pi_E OR by ' +
                '(pi_E_N and pi_E_E).')

        # If parallax is defined, then both components must be set:
        if ('pi_E_N' in keys) != ('pi_E_E' in keys):
            raise KeyError(
                'You have to define either both or none of (pi_E_N, pi_E_E).')

        # t_0_par makes sense only when parallax is defined.
        if 't_0_par' in keys and not self._type['parallax']:
            raise KeyError('t_0_par makes sense only when parallax is defined')

        # Parallax needs reference epoch:
        if 'pi_E' in keys or 'pi_E_N' in keys:
            if 't_0' not in keys and 't_0_par' not in keys:
                raise KeyError(
                    'Parallax is defined, hence either t_0 or t_0_par has ' +
                    'to be set.')

    def _check_valid_combination_1_source_mass_sheet(self, keys):
        """
        Make sure that alpha is defined if shear_G is defined,
        but not if only convergence_K is defined.
        """
        if ('shear_G' in keys) and ('alpha' not in keys):
            raise KeyError(
                'A model with external mass sheet shear requires alpha.')

        if ('shear_G' not in keys) and ('convergence_K' in keys):
            if 'alpha' in keys:
                raise KeyError(
                    'A model with external mass sheet convergence and '
                    'no shear cannot have alpha defined.')

    def _check_valid_combination_1_source_binary_lens(self, keys):
        """
        Here we check binary lens parameters for non-Cassan08 parameterization.
        """
        # s, q, and alpha must all be defined if s or q are defined
        if ('s' in keys) or ('q' in keys):
            if (('s' not in keys) or
                    ('q' not in keys) or ('alpha' not in keys)):
                raise KeyError(
                    'A binary model requires all three of (s, q, alpha).')

        # If ds_dt is defined, dalpha_dt must be defined
        if ('ds_dt' in keys) or ('dalpha_dt' in keys):
            if ('ds_dt' not in keys) or ('dalpha_dt' not in keys):
                raise KeyError(
                    'Lens orbital motion requires both ds_dt and dalpha_dt.' +
                    '\nNote that you can set either of them to 0.')
        # If orbital motion is defined, then reference epoch has to be set.
            if 't_0' not in keys and 't_0_kep' not in keys:
                raise KeyError(
                    'Orbital motion requires reference epoch, ' +
                    'i.e., t_0 or t_0_kep')

        # t_0_kep makes sense only when orbital motion is defined.
        if 't_0_kep' in keys:
            if 'ds_dt' not in keys or 'dalpha_dt' not in keys:
                raise KeyError(
                    't_0_kep makes sense only when orbital motion is defined.')

    def _check_valid_combination_1_source_xallarap(self, keys):
        """
        If xallarap parameters are defined,
        then make sure there are all required parameters
        """
        if not self._type['xallarap']:
            return

        self._check_orbit_parameters(keys, "xi_")

    def _check_orbit_parameters(self, keys, prefix):
        """
        check if orbit is properly defined; prefix is added to
        checked orbit parameters
        """
        required = ('period semimajor_axis inclination '
                    'Omega_node argument_of_latitude_reference').split()
        required = [prefix + req for req in required]
        for parameter in required:
            if parameter not in keys:
                raise KeyError(parameter)

        allowed = set([prefix + 'eccentricity', prefix + 'omega_periapsis'])
        n_used = len(set(keys).intersection(allowed))
        if n_used not in [0, len(allowed)]:
            raise KeyError(
                'Error in defining ' + prefix + 'eccentricity and ' +
                prefix + 'omega_periapsis. ' +
                'Both of them or neither should be defined.')

    def _check_valid_combination_1_source_Cassan08(self, keys):
        """
        Check parameters defined for Cassan 2008 parameterization.
        Currently, only static models are accepted.
        """
        # Check that all required parameters are defined.
        parameters = ['s', 'q', 'x_caustic_in', 'x_caustic_out',
                      't_caustic_in', 't_caustic_out']
        for parameter in parameters:
            if parameter not in keys:
                raise KeyError(
                    'If you use Cassan 2008 parameterization, then all ' +
                    'these parameters have to be defined:\n' +
                    ' \n'.join(parameters))

        # Source size cannot be over-defined.
        if ('rho' in keys) and ('t_star' in keys):
            raise KeyError('Both rho and t_star cannot be defined for ' +
                           'Cassan 08 parameterization.')

    def _check_valid_parameter_values(self, parameters):
        """
        Prevent user from setting negative (unphysical) values for
        t_E, t_star, rho etc. Shear_G should be complex.

        Also, check that all values are scalars (except pi_E vector).
        """
        full_names = {
            't_E': 'Einstein timescale', 't_star': 'Source crossing time',
            'rho': 'Source size', 's': 'separation'}

        for (name, full) in full_names.items():
            if name in parameters.keys():
                if parameters[name] < 0.:
                    fmt = "{:} cannot be negative: {:}"
                    raise ValueError(fmt.format(full, parameters[name]))

        for (key, value) in parameters.items():
            if key == 'pi_E':
                continue
            check = (not np.isscalar(value) or isinstance(value, str))
            if not isinstance(value, u.Quantity) and check:
                msg = "{:} must be a scalar: {:}, {:}"
                raise TypeError(msg.format(key, value, type(value)))

        for name in ['x_caustic_in', 'x_caustic_out']:
            if name in parameters.keys():
                if parameters[name] < 0. or parameters[name] > 1.:
                    msg = "Parameter {:} has to be in [0, 1] range, not {:}"
                    raise ValueError(msg.format(name, parameters[name]))

        for name in ['q']:
            if name in parameters.keys():
                if parameters[name] <= 0. or parameters[name] >= 1.:
                    msg = "Parameter {:} has to be in (0, 1) range, not {:}"
                    raise ValueError(msg.format(name, parameters[name]))

        for name in ['xi_eccentricity']:
            if name in parameters.keys():
                if parameters[name] < 0. or parameters[name] >= 1.:
                    msg = "Parameter {:} has to be in [0, 1) range, not {:}"
                    raise ValueError(msg.format(name, parameters[name]))

        if 'shear_G' in parameters.keys():
            if not isinstance(parameters['shear_G'], complex):
                raise TypeError("External shear (shear_G) must be complex")

    def _set_parameters(self, parameters):
        """
        check if parameter values make sense and remember the copy of the dict
        """
        self._check_valid_parameter_values(parameters)
        self.parameters = dict(parameters)

        for parameter in ['t_E', 't_star', 't_eff', 't_star_1', 't_star_2']:
            if parameter in self.parameters:
                self._set_time_quantity(parameter, self.parameters[parameter])

        angle_parameters = [
            'alpha', 'xi_Omega_node', 'xi_inclination',
            'xi_argument_of_latitude_reference', 'xi_omega_periapsis']
        for parameter in angle_parameters:
            if parameter in self.parameters:
                self._warn_if_angle_outside_reasonable_range(
                    self.parameters[parameter], parameter)

    def _update_sources(self, parameter, value):
        """
        For multi-source models, update the values for all sources.
        Note that pi_E_N and pi_E_E are changed separately.
        """
        if self.n_sources == 1:
            return

        if parameter in self._source_1_parameters.parameters:
            setattr(self._source_1_parameters, parameter, value)

        if parameter in self._source_2_parameters.parameters:
            setattr(self._source_2_parameters, parameter, value)

    def _set_time_quantity(self, key, new_time):
        """
        Save a variable with units of time (e.g. t_E, t_star,
        t_eff). If units are not given, assume days.
        """
        if isinstance(new_time, u.Quantity):
            self.parameters[key] = new_time
        else:
            self.parameters[key] = new_time * u.day

    def _check_time_quantity(self, key):
        """
        Make sure that value for give key has quantity, add it if missing.
        """
        if not isinstance(self.parameters[key], u.Quantity):
            self._set_time_quantity(key, self.parameters[key])

    def _get_uniform_caustic_sampling(self):
        """
        Sets self._uniform_caustic if that is required.
        Also resets self._standard_parameters.
        """
        recalculate = (self._uniform_caustic is None or
                       self.s != self._uniform_caustic.s or
                       self.q != self._uniform_caustic.q)
        if recalculate:
            self._uniform_caustic = UniformCausticSampling(s=self.s, q=self.q)
            self._standard_parameters = None

    def _get_standard_parameters_from_Cassan08(self):
        """
        Calculate these parameters:
        t_0 u_0 t_E alpha
        based on:
        x_caustic_in x_caustic_out t_caustic_in t_caustic_out
        using transformation that depends on:
        s q
        """
        self._get_uniform_caustic_sampling()

        if self._standard_parameters is None:
            keys = ['x_caustic_in', 'x_caustic_out',
                    't_caustic_in', 't_caustic_out']
            kwargs = {key: self.parameters[key] for key in keys}
            self._standard_parameters = (
                self._uniform_caustic.get_standard_parameters(**kwargs))

    @property
    def t_0(self):
        """
        *float*

        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['t_0']
        return self.parameters['t_0']

    @t_0.setter
    def t_0(self, new_t_0):
        if self._type['Cassan08']:
            raise ValueError('t_0 cannot be set for model using ' +
                             'Cassan (2008) parameterization')
        self.parameters['t_0'] = new_t_0
        self._update_sources('t_0', new_t_0)

    @property
    def u_0(self):
        """
        *float*

        The minimum projected separation between the source
        and the lens center of mass.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['u_0']
        if 'u_0' in self.parameters.keys():
            return self.parameters['u_0']
        else:
            try:
                u_0_quantity = (
                    self.parameters['t_eff'] / self.parameters['t_E'])
                return (u_0_quantity + 0.).value
                # Adding 0 ensures the units are simplified.
            except KeyError:
                raise AttributeError(
                    'u_0 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @u_0.setter
    def u_0(self, new_u_0):
        if self._type['Cassan08']:
            raise ValueError('u_0 cannot be set for model using ' +
                             'Cassan (2008) parameterization')
        if 'u_0' in self.parameters.keys():
            self.parameters['u_0'] = new_u_0
            self._update_sources('u_0', new_u_0)
        else:
            raise KeyError('u_0 is not a parameter of this model.')

    @property
    def t_star(self):
        """
        *float*

        t_star = rho * t_E = source radius crossing time

        "day" is the default unit. Can be set as *float* or
        *astropy.Quantity*, but always returns *float* in units of days.
        """
        if 't_star' in self.parameters.keys():
            self._check_time_quantity('t_star')
            return self.parameters['t_star'].to(u.day).value
        elif ('rho' in self.parameters.keys() and self._type['Cassan08']):
            return self.rho * self.t_E
        else:
            try:
                return (self.parameters['t_E'].to(u.day).value *
                        self.parameters['rho'])
            except KeyError:
                raise AttributeError(
                    't_star is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star.setter
    def t_star(self, new_t_star):
        if 't_star' in self.parameters.keys():
            self._set_time_quantity('t_star', new_t_star)
            self._update_sources('t_star', new_t_star)
        else:
            raise KeyError('t_star is not a parameter of this model.')

        if new_t_star < 0.:
            raise ValueError(
                'Source crossing time cannot be negative:', new_t_star)

    @property
    def t_eff(self):
        """
        *float*

        t_eff = u_0 * t_E = effective timescale

        "day" is the default unit. Can be set as *float* or
        *astropy.Quantity*, but always returns *float* in units of days.
        """
        if 't_eff' in self.parameters.keys():
            self._check_time_quantity('t_eff')
            return self.parameters['t_eff'].to(u.day).value
        else:
            try:
                return (self.parameters['t_E'].to(u.day).value *
                        self.parameters['u_0'])
            except KeyError:
                raise AttributeError(
                    't_eff is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_eff.setter
    def t_eff(self, new_t_eff):
        if 't_eff' in self.parameters.keys():
            self._set_time_quantity('t_eff', new_t_eff)
            self._update_sources('t_eff', new_t_eff)
        else:
            raise KeyError('t_eff is not a parameter of this model.')

    @property
    def t_E(self):
        """
        *float*

        The Einstein timescale. "day" is the default unit. Can be set as
        *float* or *astropy.Quantity*, but always returns *float* in units of
        days.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['t_E']
        if 't_E' in self.parameters.keys():
            self._check_time_quantity('t_E')
            return self.parameters['t_E'].to(u.day).value
        elif ('t_star' in self.parameters.keys() and
              'rho' in self.parameters.keys()):
            return self.t_star / self.rho
        elif ('t_eff' in self.parameters.keys() and
              'u_0' in self.parameters.keys()):
            return self.t_eff / abs(self.u_0)
        else:
            raise KeyError("You're trying to access t_E that was not set")

    @t_E.setter
    def t_E(self, new_t_E):
        if self._type['Cassan08']:
            raise ValueError('t_E cannot be set for model using ' +
                             'Cassan (2008) parameterization')

        if new_t_E is None:
            raise ValueError('Must provide a value')

        if new_t_E < 0.:
            raise ValueError('Einstein timescale cannot be negative:', new_t_E)

        if 't_E' in self.parameters.keys():
            self._set_time_quantity('t_E', new_t_E)
            self._update_sources('t_E', new_t_E)
        else:
            raise KeyError('t_E is not a parameter of this model.')

    @property
    def rho(self):
        """
        *float*

        source size as a fraction of the Einstein radius
        """
        if 'rho' in self.parameters.keys():
            return self.parameters['rho']
        elif ('t_star' in self.parameters.keys() and
              't_E' in self.parameters.keys()):
            return self.t_star / self.t_E
        elif ('t_star' in self.parameters.keys() and self._type['Cassan08']):
            return self.t_star / self.t_E
        else:
            return None

    @rho.setter
    def rho(self, new_rho):
        if 'rho' in self.parameters.keys():
            if new_rho < 0.:
                raise ValueError('source size (rho) cannot be negative')
            self.parameters['rho'] = new_rho
            self._update_sources('rho', new_rho)
        else:
            raise KeyError('rho is not a parameter of this model.')

    @property
    def alpha(self):
        """
        *astropy.Quantity*

        The angle of the source trajectory relative to the binary lens
        axis (or primary-secondary axis). Measured counterclockwise,
        i.e., according to convention advocated by
        `Skowron et al. 2011 (ApJ, 738, 87)
        <https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract>`_,
        but shifted by 180 deg.  May be
        set as a *float* --> assumes "deg" is the default unit.
        Regardless of input value, returns value in degrees.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['alpha'] * u.deg

        if not isinstance(self.parameters['alpha'], u.Quantity):
            self.parameters['alpha'] = self.parameters['alpha'] * u.deg
        return self.parameters['alpha'].to(u.deg)

    @alpha.setter
    def alpha(self, new_alpha):
        if self._type['Cassan08']:
            raise ValueError('alpha cannot be set for model using ' +
                             'Cassan (2008) parameterization')

        if isinstance(new_alpha, u.Quantity):
            self.parameters['alpha'] = new_alpha
        else:
            self.parameters['alpha'] = new_alpha * u.deg
        self._warn_if_angle_outside_reasonable_range(
            self.parameters['alpha'].to(u.deg).value, 'alpha')
        self._update_sources('alpha', new_alpha)

    def _warn_if_angle_outside_reasonable_range(self, value, name):
        """
        Check if value of given angle is in reasonable range and warn if not
        """
        min_ = -360.
        max_ = 540.
        if isinstance(value, u.Quantity):
            value = value.to(u.deg).value
        if value < min_ or value > max_:
            fmt = "Strange value of angle {:}: {:}"
            warnings.warn(fmt.format(name, value), RuntimeWarning)

    @property
    def q(self):
        """
        *float*

        mass ratio of the two lens components. Only 2 bodies allowed.
        """
        return self.parameters['q']

    @q.setter
    def q(self, new_q):
        if new_q < 0. or new_q > 1.:
            raise ValueError('mass ratio q has to be between 0 and 1')
        self.parameters['q'] = new_q
        self._update_sources('q', new_q)

    @property
    def convergence_K(self):
        """
        *float*

        Convergence of external mass sheet.
        """
        if 'convergence_K' in self.parameters.keys():
            return self.parameters['convergence_K']
        else:
            raise KeyError('convergence_K is not a parameter of this model.')

    @convergence_K.setter
    def convergence_K(self, new_K):
        if 'convergence_K' in self.parameters.keys():
            self.parameters['convergence_K'] = new_K
            self._update_sources('convergence_K', new_K)
        else:
            raise KeyError('convergence_K is not a parameter of this model.')

    @property
    def shear_G(self):
        """
        *complex*

        Shear of external mass sheet.
        """
        if 'shear_G' in self.parameters.keys():
            return self.parameters['shear_G']
        else:
            return None

    @shear_G.setter
    def shear_G(self, new_G):
        if 'shear_G' in self.parameters.keys():
            self.parameters['shear_G'] = new_G
            self._update_sources('shear_G', new_G)
        else:
            raise KeyError('shear_G is not a parameter of this model.')

    @property
    def s(self):
        """
        *float*

        separation of the two lens components relative to Einstein ring size
        """
        return self.parameters['s']

    @s.setter
    def s(self, new_s):
        if new_s < 0.:
            raise ValueError(
                'Binary lens separation cannot be negative:', new_s)

        self.parameters['s'] = new_s
        self._update_sources('s', new_s)

    @property
    def pi_E(self):
        """
        *list of floats*

        The microlensing parallax vector. Must be set as a vector/list
        (i.e. [pi_E_N, pi_E_E]). To get the magnitude of pi_E, use
        pi_E_mag
        """
        if 'pi_E' in self.parameters.keys():
            return self.parameters['pi_E']
        elif ('pi_E_N' in self.parameters.keys() and
              'pi_E_E' in self.parameters.keys()):
            return [self.parameters['pi_E_N'], self.parameters['pi_E_E']]
        else:
            return None

    @pi_E.setter
    def pi_E(self, new_pi_E):
        if isinstance(new_pi_E, np.ndarray):
            new_pi_E = new_pi_E.flatten()

        if 'pi_E' in self.parameters.keys():
            if len(new_pi_E) == 2:
                self.parameters['pi_E'] = new_pi_E
                self._update_sources('pi_E', new_pi_E)
            else:
                raise TypeError('pi_E is a 2D vector. It must have length 2.')

        elif ('pi_E_N' in self.parameters.keys() and
              'pi_E_E' in self.parameters.keys()):
            self.parameters['pi_E_N'] = new_pi_E[0]
            self.parameters['pi_E_E'] = new_pi_E[1]
            self._update_sources('pi_E_N', new_pi_E[0])
            self._update_sources('pi_E_E', new_pi_E[1])
        else:
            raise KeyError('pi_E is not a parameter of this model.')

    @property
    def pi_E_N(self):
        """
        *float*

        The North component of the microlensing parallax vector.
        """
        if 'pi_E_N' in self.parameters.keys():
            return self.parameters['pi_E_N']
        elif 'pi_E' in self.parameters.keys():
            return self.parameters['pi_E'][0]
        else:
            raise KeyError('pi_E_N not defined for this model')

    @pi_E_N.setter
    def pi_E_N(self, new_value):
        if 'pi_E_N' in self.parameters.keys():
            self.parameters['pi_E_N'] = new_value
            self._update_sources('pi_E_N', new_value)
        elif 'pi_E' in self.parameters.keys():
            self.parameters['pi_E'][0] = new_value
            if self.n_sources != 1:
                self._source_1_parameters.parameters['pi_E'][0] = new_value
                self._source_2_parameters.parameters['pi_E'][0] = new_value
        else:
            raise KeyError('pi_E_N is not a parameter of this model.')

    @property
    def pi_E_E(self):
        """
        *float*

        The East component of the microlensing parallax vector.
        """
        if 'pi_E_E' in self.parameters.keys():
            return self.parameters['pi_E_E']
        elif 'pi_E' in self.parameters.keys():
            return self.parameters['pi_E'][1]
        else:
            raise KeyError('pi_E_N not defined for this model')

    @pi_E_E.setter
    def pi_E_E(self, new_value):
        if 'pi_E_E' in self.parameters.keys():
            self.parameters['pi_E_E'] = new_value
            self._update_sources('pi_E_E', new_value)
        elif 'pi_E' in self.parameters.keys():
            self.parameters['pi_E'][1] = new_value
            if self.n_sources != 1:
                self._source_1_parameters.parameters['pi_E'][1] = new_value
                self._source_2_parameters.parameters['pi_E'][1] = new_value
        else:
            raise KeyError('pi_E_E is not a parameter of this model.')

    @property
    def t_0_par(self):
        """
        *float*

        The reference time for the calculation of parallax. If not set
        explicitly, then it is assumed t_0_par = t_0.

        Note that this is a reference value and not the fitting parameter.
        It is best to fix it at the begin of calculations.
        """
        if 't_0_par' not in self.parameters.keys():
            return self.parameters['t_0']
        else:
            return self.parameters['t_0_par']

    @t_0_par.setter
    def t_0_par(self, new_t_0_par):
        self.parameters['t_0_par'] = new_t_0_par
        self._update_sources('t_0_par', new_t_0_par)

    @property
    def pi_E_mag(self):
        """
        *float*

        The magnitude of the microlensing parallax vector.
        """
        if 'pi_E' in self.parameters.keys():
            pi_E_N = self.parameters['pi_E'][0]
            pi_E_E = self.parameters['pi_E'][1]
        elif ('pi_E_N' in self.parameters.keys() and
              'pi_E_E' in self.parameters.keys()):
            pi_E_N = self.parameters['pi_E_N']
            pi_E_E = self.parameters['pi_E_E']
        else:
            raise KeyError('pi_E not defined for this model')
        return np.sqrt(pi_E_N**2 + pi_E_E**2)

    @property
    def x_caustic_in(self):
        """
        *float*

        Curvelinear coordinate (in `Cassan (2008) parameterization
        <https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract>`_)
        of caustic entrance for a static binary lens model. See
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`.
        """
        return self.parameters['x_caustic_in']

    @x_caustic_in.setter
    def x_caustic_in(self, new_value):
        if new_value < 0. or new_value > 1.:
            msg = "x_caustic_in must be between 0 and 1, not {:}"
            raise ValueError(msg.format(new_value))
        self._standard_parameters = None
        self.parameters['x_caustic_in'] = new_value

    @property
    def x_caustic_out(self):
        """
        *float*

        Curvelinear coordinate (in `Cassan (2008) parameterization
        <https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract>`_)
        of caustic exit for a static binary lens model. See
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`.
        """
        return self.parameters['x_caustic_out']

    @x_caustic_out.setter
    def x_caustic_out(self, new_value):
        if new_value < 0. or new_value > 1.:
            msg = "x_caustic_out must be between 0 and 1, not {:}"
            raise ValueError(msg.format(new_value))
        self._standard_parameters = None
        self.parameters['x_caustic_out'] = new_value

    @property
    def t_caustic_in(self):
        """
        *float*

        Epoch of caustic entrance for a static binary lens model in
        `Cassan (2008) parameterization
        <https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract>`_)
        See
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`.
        """
        return self.parameters['t_caustic_in']

    @t_caustic_in.setter
    def t_caustic_in(self, new_value):
        self._standard_parameters = None
        self.parameters['t_caustic_in'] = new_value

    @property
    def t_caustic_out(self):
        """
        *float*

        Epoch of caustic exit for a static binary lens model in
        `Cassan (2008) parameterization
        <https://ui.adsabs.harvard.edu/abs/2008A%26A...491..587C/abstract>`_)
        See
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`.
        """
        return self.parameters['t_caustic_out']

    @t_caustic_out.setter
    def t_caustic_out(self, new_value):
        self._standard_parameters = None
        self.parameters['t_caustic_out'] = new_value

    @property
    def ds_dt(self):
        """
        *astropy.Quantity*

        Change rate of separation :py:attr:`~s` in 1/year. Can be set as
        *AstroPy.Quantity* or as *float* (1/year is assumed default unit).
        Regardless of input value, returns value in 1/year.
        """
        if not isinstance(self.parameters['ds_dt'], u.Quantity):
            self.parameters['ds_dt'] = self.parameters['ds_dt'] / u.yr

        return self.parameters['ds_dt'].to(1 / u.yr)

    @ds_dt.setter
    def ds_dt(self, new_ds_dt):
        if isinstance(new_ds_dt, u.Quantity):
            self.parameters['ds_dt'] = new_ds_dt
        else:
            self.parameters['ds_dt'] = new_ds_dt / u.yr
        self._update_sources('ds_dt', new_ds_dt)

    @property
    def dalpha_dt(self):
        """
        *astropy.Quantity*

        Change rate of angle :py:attr:`~alpha` in deg/year. Can be set as
        *AstroPy.Quantity* or as *float* (deg/year is assumed default unit).
        Regardless of input value, returns value in deg/year.
        """
        if not isinstance(self.parameters['dalpha_dt'], u.Quantity):
            self.parameters['dalpha_dt'] = (self.parameters['dalpha_dt'] *
                                            u.deg / u.yr)

        return self.parameters['dalpha_dt'].to(u.deg / u.yr)

    @dalpha_dt.setter
    def dalpha_dt(self, new_dalpha_dt):
        if isinstance(new_dalpha_dt, u.Quantity):
            self.parameters['dalpha_dt'] = new_dalpha_dt
        else:
            self.parameters['dalpha_dt'] = new_dalpha_dt * u.deg / u.yr
        self._update_sources('dalpha_dt', new_dalpha_dt)

    @property
    def t_0_kep(self):
        """
        *float*

        The reference time for the calculation of lens orbital motion.
        If not set explicitly, then it is assumed t_0_kep = t_0.

        Note that this is a reference value and not the fitting parameter.
        It is best to fix it at the begin of calculations.
        """
        if 't_0_kep' not in self.parameters.keys():
            return self.parameters['t_0']
        else:
            return self.parameters['t_0_kep']

    @t_0_kep.setter
    def t_0_kep(self, new):
        self.parameters['t_0_kep'] = new
        self._update_sources('t_0_kep', new)

    @property
    def xi_period(self):
        """
        *float*

        Orbital period of the source system (xallarap) in days.
        """
        return self.parameters['xi_period']

    @xi_period.setter
    def xi_period(self, new_value):
        if new_value < 0.:
            raise ValueError('Xallarap period cannot be negative')
        self.parameters['xi_period'] = new_value

    @property
    def xi_semimajor_axis(self):
        """
        *float*

        Semi-major axis of the source orbit (xallarap) in the theta_E units.
        """
        return self.parameters['xi_semimajor_axis']

    @xi_semimajor_axis.setter
    def xi_semimajor_axis(self, new_value):
        if new_value < 0.:
            raise ValueError('Xallarap semimajor axis cannot be negative')
        self.parameters['xi_semimajor_axis'] = new_value

    @property
    def xi_Omega_node(self):
        """
        *float*

        The longitude of the ascending node of the xallarap orbit, i.e.,
        the angle from relative lens-source proper motion direction
        to the ascending node direction.
        The units are degrees.
        """
        return self.parameters['xi_Omega_node']

    @xi_Omega_node.setter
    def xi_Omega_node(self, new_value):
        self._warn_if_angle_outside_reasonable_range(new_value,
                                                     'xi_Omega_node')
        self.parameters['xi_Omega_node'] = new_value

    @property
    def xi_inclination(self):
        """
        *float*

        The inclination of the xallarap orbit, i.e.,
        the angle between source-orbit plane and the sky plane.
        The units are degrees.
        """
        return self.parameters['xi_inclination']

    @xi_inclination.setter
    def xi_inclination(self, new_value):
        self._warn_if_angle_outside_reasonable_range(new_value,
                                                     'xi_inclination')
        self.parameters['xi_inclination'] = new_value

    @property
    def xi_argument_of_latitude_reference(self):
        """
        *float*

        The argument of latitude for the xallarap orbit at :py:attr:`~t_0_xi`.
        The argument of latitude is a sum of the true anomaly and
        the argument of periapsis. In standard notation: u = nu + omega.
        This parameter is internally used to calculate perihelion passage
        (T_0 in standard notation).
        The units are degrees.
        """
        return self.parameters['xi_argument_of_latitude_reference']

    @xi_argument_of_latitude_reference.setter
    def xi_argument_of_latitude_reference(self, new_value):
        self._warn_if_angle_outside_reasonable_range(
            new_value, 'xi_argument_of_latitude_reference')
        self.parameters['xi_argument_of_latitude_reference'] = new_value

    @property
    def xi_eccentricity(self):
        """
        *float*

        The eccentricity of the xallarap orbit. Has to be in [0, 1) range.
        """
        return self.parameters['xi_eccentricity']

    @xi_eccentricity.setter
    def xi_eccentricity(self, new_value):
        if new_value < 0. or new_value >= 1.:
            raise ValueError('xallarap eccentricity has to be between 0 and 1')
        self.parameters['xi_eccentricity'] = new_value

    @property
    def xi_omega_periapsis(self):
        """
        *float*

        The argument of periapsis of the xallrap orbit, i.e., the angle
        between the ascending node and periapsis measured in
        the direction of motion.
        The units are degrees.
        """
        return self.parameters['xi_omega_periapsis']

    @xi_omega_periapsis.setter
    def xi_omega_periapsis(self, new_value):
        self._warn_if_angle_outside_reasonable_range(
            new_value, 'xi_omega_periapsis')
        self.parameters['xi_omega_periapsis'] = new_value

    @property
    def t_0_xi(self):
        """
        *float*

        Reference epoch for xallarap orbit.
        If not provided, then it defaults to :py:attr:`~t_0`.
        """
        if 't_0_xi' not in self.parameters.keys():
            return self.parameters['t_0']
        else:
            return self.parameters['t_0_xi']

    @t_0_xi.setter
    def t_0_xi(self, new_value):
        self.parameters['t_0_xi'] = new_value

    @property
    def t_0_1(self):
        """
        *float*

        The time of minimum projected separation between the source no. 1
        and the lens center of mass.
        """
        return self.parameters['t_0_1']

    @t_0_1.setter
    def t_0_1(self, new_t_0_1):
        self.parameters['t_0_1'] = new_t_0_1
        self._source_1_parameters.t_0 = new_t_0_1

    @property
    def t_0_2(self):
        """
        *float*

        The time of minimum projected separation between the source no. 2
        and the lens center of mass.
        """
        return self.parameters['t_0_2']

    @t_0_2.setter
    def t_0_2(self, new_t_0_2):
        self.parameters['t_0_2'] = new_t_0_2
        self._source_2_parameters.t_0 = new_t_0_2

    @property
    def u_0_1(self):
        """
        *float*

        The minimum projected separation between the source no. 1
        and the lens center of mass.
        """
        if 'u_0_1' in self.parameters.keys():
            return self.parameters['u_0_1']
        else:
            try:
                t_eff = self._source_1_parameters.parameters['t_eff']
                t_E = self._source_1_parameters.parameters['t_E']
                return t_eff / t_E
            except KeyError:
                raise AttributeError(
                    'u_0_1 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @u_0_1.setter
    def u_0_1(self, new_u_0_1):
        if 'u_0_1' in self.parameters.keys():
            self.parameters['u_0_1'] = new_u_0_1
            self._source_1_parameters.u_0 = new_u_0_1
        else:
            raise KeyError('u_0_1 is not a parameter of this model.')

    @property
    def u_0_2(self):
        """
        *float*

        The minimum projected separation between the source no. 2
        and the lens center of mass.
        """
        if 'u_0_2' in self.parameters.keys():
            return self.parameters['u_0_2']
        else:
            try:
                t_eff = self._source_2_parameters.parameters['t_eff']
                t_E = self._source_2_parameters.parameters['t_E']
                return t_eff / t_E
            except KeyError:
                raise AttributeError(
                    'u_0_2 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @u_0_2.setter
    def u_0_2(self, new_u_0_2):
        if 'u_0_2' in self.parameters.keys():
            self.parameters['u_0_2'] = new_u_0_2
            self._source_2_parameters.u_0 = new_u_0_2
        else:
            raise KeyError('u_0_2 is not a parameter of this model.')

    @property
    def t_star_1(self):
        """
        *float*

        t_star_1 = rho_1 * t_E_1 = source no. 1 radius crossing time

        "day" is the default unit. Can be set as *float* or
        *astropy.Quantity*, but always returns *float* in units of days.
        """
        if 't_star_1' in self.parameters.keys():
            self._check_time_quantity('t_star_1')
            return self.parameters['t_star_1'].to(u.day).value
        else:
            try:
                t_E = self._source_1_parameters.parameters['t_E'].to(u.day)
                rho = self._source_1_parameters.parameters['rho']
                return t_E.value * rho
            except KeyError:
                raise AttributeError(
                    't_star_1 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star_1.setter
    def t_star_1(self, new_t_star_1):
        if 't_star_1' in self.parameters.keys():
            self._set_time_quantity('t_star_1', new_t_star_1)
            self._source_1_parameters.t_star = new_t_star_1
        else:
            raise KeyError('t_star_1 is not a parameter of this model.')

        if new_t_star_1 < 0.:
            raise ValueError(
                'Source crossing time cannot be negative:', new_t_star_1)

    @property
    def t_star_2(self):
        """
        *float*

        t_star_2 = rho_2 * t_E_2 = source no. 2 radius crossing time

        "day" is the default unit. Can be set as *float* or
        *astropy.Quantity*, but always returns *float* in units of days.
        """
        if 't_star_2' in self.parameters.keys():
            self._check_time_quantity('t_star_2')
            return self.parameters['t_star_2'].to(u.day).value
        else:
            try:
                t_E = self._source_2_parameters.parameters['t_E'].to(u.day)
                rho = self._source_2_parameters.parameters['rho']
                return t_E.value * rho
            except KeyError:
                raise AttributeError(
                    't_star_2 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star_2.setter
    def t_star_2(self, new_t_star_2):
        if 't_star_2' in self.parameters.keys():
            self._set_time_quantity('t_star_2', new_t_star_2)
            self._source_2_parameters.t_star = new_t_star_2
        else:
            raise KeyError('t_star_2 is not a parameter of this model.')

        if new_t_star_2 < 0.:
            raise ValueError(
                'Source crossing time cannot be negative:', new_t_star_2)

    @property
    def rho_1(self):
        """
        *float*

        source no. 1 size as a fraction of the Einstein radius
        """
        if 'rho_1' in self.parameters.keys():
            return self.parameters['rho_1']
        elif ('t_star' in self._source_1_parameters.parameters.keys() and
                't_E' in self._source_1_parameters.parameters.keys()):
            return (self._source_1_parameters.t_star /
                    self._source_1_parameters.t_E)
        else:
            return None

    @rho_1.setter
    def rho_1(self, new_rho_1):
        if 'rho_1' in self.parameters.keys():
            if new_rho_1 < 0.:
                raise ValueError('source size (rho_1) cannot be negative')
            self.parameters['rho_1'] = new_rho_1
            self._source_1_parameters.rho = new_rho_1
        else:
            raise KeyError('rho_1 is not a parameter of this model.')

    @property
    def rho_2(self):
        """
        *float*

        source no. 2 size as a fraction of the Einstein radius
        """
        if 'rho_2' in self.parameters.keys():
            return self.parameters['rho_2']
        elif ('t_star' in self._source_2_parameters.parameters.keys() and
                't_E' in self._source_2_parameters.parameters.keys()):
            return (self._source_2_parameters.t_star /
                    self._source_2_parameters.t_E)
        else:
            return None

    @rho_2.setter
    def rho_2(self, new_rho_2):
        if 'rho_2' in self.parameters.keys():
            if new_rho_2 < 0.:
                raise ValueError('source size (rho_2) cannot be negative')
            self.parameters['rho_2'] = new_rho_2
            self._source_2_parameters.rho = new_rho_2
        else:
            raise KeyError('rho_2 is not a parameter of this model.')

    def get_s(self, epoch):
        """
        Returns the value of separation :py:attr:`~s` at a given epoch or
        epochs (if orbital motion parameters are set).

        Arguments :
            epoch: *float*, *list*, *np.ndarray*
                The time(s) at which to calculate :py:attr:`~s`.

        Returns :
            separation: *float* or *np.ndarray*
                Value(s) of separation for given epochs.

        """
        if 'ds_dt' not in self.parameters.keys():
            return self.s

        if isinstance(epoch, list):
            epoch = np.array(epoch)

        s_of_t = (self.s + self.ds_dt * (epoch - self.t_0_kep) * u.d).value

        return s_of_t

    def get_alpha(self, epoch):
        """
        Returns the value of angle :py:attr:`~alpha` at a given epoch or
        epochs (if orbital motion parameters are set).

        Arguments :
            epoch: *float*, *list*, *np.ndarray*
                The time(s) at which to calculate :py:attr:`~alpha`.

        Returns :
            separation: *astropy.Quantity*
                Value(s) of angle for given epochs in degrees

        """
        if 'dalpha_dt' not in self.parameters.keys():
            return self.alpha

        if isinstance(epoch, list):
            epoch = np.array(epoch)

        alpha_of_t = (self.alpha + self.dalpha_dt *
                      (epoch - self.t_0_kep) * u.d)

        return alpha_of_t.to(u.deg)

    @property
    def gamma_parallel(self):
        """
        *astropy.Quantity*

        Parallel component of instantaneous velocity of the secondary
        relative to the primary in 1/year.
        It is parallel to the primary-secondary axis.
        Equals :py:attr:`~ds_dt`/:py:attr:`~s`. Cannot be set.
        """
        return self.ds_dt / self.s

    @property
    def gamma_perp(self):
        """
        *astropy.Quantity*

        Perpendicular component of instantaneous velocity of the secondary
        relative to the primary. It is perpendicular to the primary-secondary
        axis. It has sign opposite to :py:attr:`~dalpha_dt`
        and is in rad/yr, not deg/yr. Cannot be set.
        """
        return -self.dalpha_dt.to(u.rad / u.yr)

    @property
    def gamma(self):
        """
        *astropy.Quantity*

        Instantaneous velocity of the secondary relative to the primary in
        1/year. Cannot be set.
        """
        gamma_perp = (self.gamma_perp / u.rad).to(1 / u.yr)
        return (self.gamma_parallel**2 + gamma_perp**2)**0.5

    def is_finite_source(self):
        """
        Checks if model has finite source. For binary source models it checks
        if either of the sources is finite.

        Returns:
            is_finite_source: *boolean*
                *True* if at least one source has finite size.
        """
        return self._type['finite source']

    def is_static(self):
        """
        Checks if model is static, i.e., orbital motion parameters are not set.

        Returns :
            is_static: *boolean*
                *True* if *dalpha_dt* or *ds_dt* are set.

        """
        return not self._type['lens 2-parameter orbital motion']

    @property
    def n_lenses(self):
        """
        *int*

        Number of objects in the lens system.
        """
        return self._n_lenses

    @property
    def n_sources(self):
        """
        *int*

        Number of luminous sources.
        It can be be 1 for a xallarap model.
        """
        return self._n_sources

    @property
    def is_external_mass_sheet(self):
        """
        *bool*

        Whether an external mass sheet is included in the model
        """
        return self._type['mass sheet']

    @property
    def is_external_mass_sheet_with_shear(self):
        """
        *bool*

        Whether an external mass sheet is included in the
        model with non-zero shear
        """
        return (('shear_G' in self.parameters.keys()) and
                (self.parameters['shear_G'] != 0))

    @property
    def is_xallarap(self):
        """
        *bool*

        Whether the parameters include the xallarap or not.
        """
        return self._type['xallarap']

    @property
    def source_1_parameters(self):
        """
        :py:class:`~MulensModel.modelparameters.ModelParameters`

        Parameters of source 1 in multi-source model.

        **Do not change returned values.** To change
        parameters of the source 1, simply change the parameters of double
        source instance.
        """
        if self.n_sources == 1:
            raise ValueError('source_1_parameters cannot be accessed for ' +
                             'single source models')
        return self._source_1_parameters

    @property
    def source_2_parameters(self):
        """
        :py:class:`~MulensModel.modelparameters.ModelParameters`

        Parameters of source 2 in multi-source model.

        **Do not change returned values.** To change
        parameters of the source 1, simply change the parameters of double
        source instance.
        """
        if self.n_sources == 1:
            raise ValueError('source_2_parameters cannot be accessed for ' +
                             'single source models')
        return self._source_2_parameters

    @property
    def uniform_caustic_sampling(self):
        """
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`

        An instance of the class
        :py:class:`~MulensModel.uniformcausticsampling.UniformCausticSampling`
        that is used to calculate standard parameters based on
        the curvelinear coordinates.
        The main usage is access to the *jacobian()* function.
        In most cases, you do not need to access this property directly.
        """
        if not self._type['Cassan08']:
            raise ValueError(
                'These parameters are not in curvelinear parameterization. ' +
                'Hence you cannot access uniform_caustic_sampling property.')

        self._get_uniform_caustic_sampling()

        return self._uniform_caustic

    def as_dict(self):
        """
        Give parameters as a dict.

        Returns :
            dictionary: *dict*
                The dictionary of model parameters.
        """
        return self.parameters
