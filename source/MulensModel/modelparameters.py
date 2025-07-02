import numpy as np

from MulensModel.uniformcausticsampling import UniformCausticSampling
from MulensModel.orbits.orbit import Orbit
from MulensModel.utils import Utils


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
            ``model_parameters.u_0 = 0.1`` or ``setattr(model_parameters, 'u_0', 0.1)``.

    Example:
        Define a point lens model:
            ``params = ModelParameters({'t_0': 2450000., 'u_0': 0.3, 't_E': 35.})``

        Then you can print the parameters:
            ``print(params)``

    """

    # parameters that may be defined for a given source
    _primary_source_params_head = ['t_0', 'u_0', 't_eff']
    _finite_source_params_head = ['rho', 't_star']
    _all_source_params_head = _primary_source_params_head + _finite_source_params_head
    _t_0_ref_types = ['par', 'kep', 'xi']

    def __init__(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('ModelParameters must be initialized with dict as a parameter\ne.g., '
                            "ModelParameters({'t_0': 2456789.0, 'u_0': 0.123, 't_E': 23.45})")

        self._count_sources(parameters.keys())
        self._count_lenses(parameters.keys())
        self._set_type(parameters.keys())
        self._check_types('alpha' in parameters.keys())

        if self.n_sources == 1:
            self._check_valid_combination_1_source(parameters.keys())
            if self._type['Cassan08']:
                self._uniform_caustic = None
                self._standard_parameters = None

            if self.is_xallarap:
                self._check_for_extra_source_parameters(parameters.keys())
                delta_1 = self._get_xallarap_position(parameters)
                self._xallarap_reference_position = delta_1

        elif self.n_sources > 1:
            self._check_valid_combination_of_sources(parameters.keys())
            if 't_E' not in parameters.keys():
                raise KeyError('Currently, the binary source calculations ' +
                               'require t_E to be directly defined, i.e., ' +
                               'has to be the same for both sources.')

            source_params = self._divide_parameters(parameters)
            for (i, params_i) in enumerate(source_params):
                # This try/except block forces checks from ._init_1_source()
                # to be run on each source parameters separately.
                try:
                    self.__setattr__('_source_{0}_parameters'.format(i + 1), ModelParameters(params_i))
                except Exception:
                    print("ERROR IN INITIALIZING SOURCE {0}".format(i + 1))
                    raise

            if self.is_xallarap:
                self._update_sources_xallarap_reference()

        else:
            msg = 'wrong number of sources. Your parameters: {0:}'
            raise ValueError(msg.format(parameters))

        self._set_parameters(parameters)

    def _update_sources_xallarap_reference(self):
        """
        Update .xallarap_reference_position for each source parameters

        Note: below we're calling private function and set private
        properties NOT of self, but self._source_X_parameters,
        which both are of the same type as self.
        """
        if self.n_sources == 1:
            self._xallarap_reference_position = self._get_xallarap_position()
        elif self.n_sources == 2:
            delta_1 = self._source_1_parameters._get_xallarap_position()
            self._source_1_parameters._xallarap_reference_position = delta_1
            self._source_2_parameters._xallarap_reference_position = delta_1
        else:
            raise ValueError('internal error')

    def _get_xallarap_position(self, parameters=None):
        """
        Get position at t_0_xi from xallarap Orbit object.

        Note: this function is called in 2 different ways:
        - directly, i.e., self._get_xallarap_orbit(), and
        - indirectly, i.e., self._source_1_parameters._get_xallarap_orbit().
        """
        if parameters is None:
            parameters = self.parameters
        t_0_xi = parameters.get('t_0_xi', parameters['t_0'])

        zip_ = parameters.items()
        orbit_parameters = {key[3:]: value for (key, value) in zip_ if key[:3] == "xi_"}
        orbit_parameters['epoch_reference'] = t_0_xi
        orbit = Orbit(**orbit_parameters)
        return orbit.get_reference_plane_position([t_0_xi])

    def __getattr__(self, item):
        (head, end) = self._split_parameter_name(item)
        if end is not None:
            return self.__getattr__('_source_{:}_parameters'.format(end)).__getattribute__(head)
        elif item.startswith("source_") and item.endswith("_parameters"):
            return object.__getattribute__(self, "_" + item)
        else:
            return object.__getattribute__(self, item)

    def _split_parameter_name(self, parameter):
        """
        Split ABC_DEF_n into ABC_DEF (str) and n (int). For parameters like t_0 or rho, n is None.
        """
        end = parameter.split('_')[-1]
        if end.isnumeric() and int(end) > 0:
            head = parameter[:-len(end)-1]
            end = int(end)
        else:
            head = parameter
            end = None

        return (head, end)

    def _count_sources(self, keys):
        """
        How many luminous sources are there?
        """
        self._n_sources = 1
        for key in keys:
            n = self._split_parameter_name(key)[1]
            if n is not None:
                if n > self._n_sources:
                    self._n_sources = n

        if 'q_source' in keys:
            if self._n_sources != 1:
                if self._n_sources != 2 and ('rho_2' in keys or 't_star_2' in keys):
                    raise RuntimeError('wrong set of parametes: ' + str(keys))
            self._n_sources = 2

    def _count_lenses(self, keys):
        """How many lenses are there?"""
        self._n_lenses = 1
        if 's' in keys or 'q' in keys:
            self._n_lenses = 2
        # Both standard and Cassan08 parameterizations require s and q

    def _set_type(self, keys):
        """
        sets self._type property, which indicates what type of a model we have
        """
        types = ['finite source', 'parallax', 'Cassan08',
                 'lens 2-parameter orbital motion', 'full keplerian motion',
                 'mass sheet', 'xallarap']
        out = {type_: False for type_ in types}

        temp = {
            'finite source': 'rho t_star rho_1 rho_2 t_star_1 t_star_2',
            'parallax': 'pi_E_N pi_E_E',
            'Cassan08': 'x_caustic_in x_caustic_out t_caustic_in t_caustic_out',
            'lens 2-parameter orbital motion': 'dalpha_dt ds_dt',
            'full keplerian motion': 's_z ds_z_dt',
            'mass sheet': 'convergence_K shear_G',
            'xallarap': ('xi_period xi_semimajor_axis xi_inclination '
                         'xi_Omega_node xi_argument_of_latitude_reference '
                         'xi_eccentricity xi_omega_periapsis q_source')}

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

        # Lens orbital motion requires binary lens:
        if self._type['lens 2-parameter orbital motion'] and n_lenses == 1:
            raise KeyError('Orbital motion of the lens requires two lens '
                           'components but only one was provided.')
        # Full Keplerian motion requires binary lens:
        if self._type['full keplerian motion'] and n_lenses == 1:
            raise KeyError('Full Keplerian motion of the lens requires two '
                           'lens components but only one was provided.')

        self._check_types_for_Cassan08()

        if alpha_defined:
            if self._n_lenses == 1 and not self._type['mass sheet']:
                raise KeyError(
                    'You defined alpha for single lens model '
                    'without external mass sheet. This is not allowed.')

        if self._type['full keplerian motion']:
            self._lens_keplerian_last_input = None
            self._lens_keplerian = dict()

    def _check_valid_combination_1_source(self, keys):
        """
        Check that the user hasn't over-defined the ModelParameters.
        """
        # Make sure that there are no unwanted keys
        allowed_keys = set((
            't_0 u_0 t_E t_eff rho t_star pi_E_N pi_E_E t_0_par '
            's q alpha dalpha_dt ds_dt s_z ds_z_dt t_0_kep convergence_K shear_G '
            't_0_1 t_0_2 u_0_1 u_0_2 rho_1 rho_2 t_star_1 t_star_2 '
            'x_caustic_in x_caustic_out t_caustic_in t_caustic_out '
            'xi_period xi_semimajor_axis xi_inclination xi_Omega_node '
            'xi_argument_of_latitude_reference xi_eccentricity '
            'xi_omega_periapsis t_0_xi q_source').split())
        difference = set(keys) - allowed_keys
        if len(difference) > 0:
            derived_1 = ['gamma', 'gamma_perp', 'gamma_parallel', 'gamma_z']
            if set(keys).intersection(derived_1):
                msg = ('You cannot set gamma, gamma_perp, gamma_parallel' +
                       'or gamma_z. These are derived parameters. You can' +
                       ' set ds_dt, dalpha_dt, s_z and ds_z_dt instead.\n')
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

        types = ['parallax', 'xallarap', 'lens 2-parameter orbital motion',
                 'full keplerian motion']
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
        Divide an input dict into each source separately.
        Some of the parameters are copied to both dicts.
        """
        skipped_parameters = ['q_source']

        source_parameters = []
        for i in range(self.n_sources):
            params_i = dict()
            for (key, value) in parameters.items():
                if key in skipped_parameters:
                    continue

                (head, end) = self._split_parameter_name(key)
                if end is None:
                    params_i[key] = value
                elif end == i + 1:
                    params_i[head] = value

            source_parameters.append(params_i)

        if self.n_sources == 2 and self._type['xallarap']:
            self._set_changed_parameters_2nd_source(parameters['q_source'], source_parameters[1])

        return source_parameters

    def _set_changed_parameters_2nd_source(self, q_source, parameters_2):
        """
        For xallarap model with 2 sources, the orbit of the second source
        must have 2 parameters changed.
        Functions starts with tests of input
        """
        if q_source <= 0.:
            raise ValueError('q_source cannot be negative')

        check_keys = ['xi_semimajor_axis', 'xi_argument_of_latitude_reference']
        for key in check_keys:
            if key not in parameters_2:
                raise KeyError('xallarap model with 2 sources requires ' + key)

        parameters_2['xi_semimajor_axis'] /= q_source
        parameters_2['xi_argument_of_latitude_reference'] += 180.
        if parameters_2['xi_argument_of_latitude_reference'] > 360.:
            parameters_2['xi_argument_of_latitude_reference'] -= 360.

    def __repr__(self):
        """A nice way to represent a ModelParameters object as a string"""
        out = self._get_main_parameters_to_print()

        if self.is_xallarap:
            fmt = "\nxallarap reference position: ({:.4f}, {:.4f})"
            if self.n_sources == 1:
                source = self
            else:
                source = self._source_1_parameters
            position = source.xallarap_reference_position
            out += fmt.format(position[0, 0], position[1, 0])

        return out

    def _get_main_parameters_to_print(self):
        """
        prepare all the standard parameters to be printed
        """
        keys = self._get_keys_for_repr()
        formats = self._get_formats_dict_for_repr()
        ordered_keys = self._get_ordered_keys_for_repr()

        variables = [''] * (self._n_sources + 1)
        values = [''] * (self._n_sources + 1)
        for key in ordered_keys:
            if key not in keys:
                continue
            index = self._split_parameter_name(key)[1]
            index = 0 if index is None else index
            (full_name, value) = self._get_values_for_repr(formats[key], key)
            (fmt_1, fmt_2) = self._get_formats_for_repr(formats[key], full_name)
            variables[index] += fmt_1.format(full_name)
            values[index] += fmt_2.format(value)

        print_msg = ''
        for (i, variable) in enumerate(variables):
            if variable and values[i]:
                print_msg += "{:}\n{:}".format(variable, values[i])
            if i < self.n_sources and variables[i+1]:
                print_msg += "\n"

        return print_msg

    def _get_keys_for_repr(self):
        """
        get all the keys that will be printed
        """
        keys = set(self.parameters.keys())

        if 'pi_E_E' in keys or 'pi_E_N' in keys:
            keys |= {'t_0_par'}

        orbital_motion_keys = ('ds_dt', 'dalpha_dt', 's_z', 'ds_z_dt')
        if any(key in keys for key in orbital_motion_keys):
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
            's_z': {'width': 11, 'precision': 5},
            'ds_z_dt': {'width': 18, 'precision': 5, 'unit': '/yr',
                        'name': 'ds_z/dt'},
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
            'q_source': {'width': 12, 'precision': 8},
            't_0_xi': {'width': 13, 'precision': 5, 'unit': 'HJD'},
        }
        # Add multiple source parameters with the same settings.
        if self.n_sources > 1:
            for i in range(self.n_sources):
                for param_head in self._all_source_params_head:
                    form = formats[param_head]
                    key = '{0}_{1}'.format(param_head, i+1)
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
        basic_keys = ['t_0', 'u_0', 't_E', 'rho', 't_star']
        additional_keys = [
            'pi_E_N', 'pi_E_E', 't_0_par', 's', 'q', 'alpha',
            'convergence_K', 'shear_G', 'ds_dt', 'dalpha_dt', 's_z',
            'ds_z_dt', 't_0_kep',
            'x_caustic_in', 'x_caustic_out', 't_caustic_in', 't_caustic_out',
            'xi_period', 'xi_semimajor_axis', 'xi_inclination',
            'xi_Omega_node', 'xi_argument_of_latitude_reference',
            'xi_eccentricity', 'xi_omega_periapsis', 'q_source', 't_0_xi'
            ]

        ordered_keys = []
        if self.n_sources > 1:
            for param_head in basic_keys:
                if param_head == 't_E':
                    ordered_keys.append(param_head)
                else:
                    for i in range(self.n_sources):
                        ordered_keys.append('{0}_{1}'.format(param_head, i + 1))

        else:
            ordered_keys = basic_keys

        for key in additional_keys:
            ordered_keys.append(key)

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
        return (full_name, value)

    def _get_formats_for_repr(self, form, full_name):
        """
        Extract formats to be used by __repr__().
        """
        fmt_1 = '{:>' + str(max([form['width'], len(full_name)]))
        fmt_2 = fmt_1 + '.' + str(form['precision']) + 'f} '
        fmt_1 += '} '
        return (fmt_1, fmt_2)

    def _check_valid_combination_of_sources(self, keys):
        """
        make sure that there is no conflict between t_0 and t_0_1 etc.
        Also make sure that t_0 and u_0 are defined for all sources.
        """
        self._check_for_incompatible_source_parameters(keys)
        self._check_for_missing_source_parameters(keys)
        self._check_for_extra_source_parameters(keys)

    def _check_for_incompatible_source_parameters(self, keys):
        """
        make sure that there is no conflict between t_0 and t_0_1 etc.
        """
        for parameter in self._primary_source_params_head:
            if parameter in keys:
                # conflict between t_0 and t_0_1
                for i in range(self.n_sources):
                    if '{0}_{1}'.format(parameter, i+1) in keys:
                        msg = 'You cannot set both {:} and {:}'
                        raise KeyError(msg.format(parameter, '{0}_{1}'.format(parameter, i+1)))

        for parameter in self._finite_source_params_head:
            if parameter in keys:
                raise KeyError('You must specify which source {0} goes with'.format(parameter))

        if self.is_xallarap:
            self._check_for_parameters_incompatible_with_xallarap(keys)

    def _check_for_parameters_incompatible_with_xallarap(self, keys):
        """
        Check for additional source parameters with xallarap (bad).
        """
        for parameter in keys:
            try:
                key_num = parameter.split('_')[-1]
                if ((int(key_num) > 1) and (parameter[0:-len(key_num) - 1] in self._primary_source_params_head)):
                    msg = 'xallarap parameters cannot be mixed with {:}'
                    raise NotImplementedError(msg.format(parameter))

            except ValueError:
                pass

    def _check_for_missing_source_parameters(self, keys):
        """
        Also make sure that t_0 and u_0 are defined for all sources.
        """
        if (self.n_sources > 1) and ('q_source' not in keys):
            for i in range(self.n_sources):
                if 't_0_{0}'.format(i + 1) not in keys:
                    raise KeyError(
                        't_0_{0} is missing from parameters.'.format(i+1) +
                        'Your parameters: {0}'.format(keys))

            for i in range(self.n_sources):
                if ('u_0_{0}'.format(i + 1) not in keys) and ('t_eff_{0}'.format(i + 1) not in keys):
                    raise KeyError(
                        'Either u_0_{0} or t_eff_{0} must be specified.'.format(i+1) +
                        'Your parameters: {0}'.format(keys))

    def _check_for_extra_source_parameters(self, keys):
        """
        Check if parameters have been set for sources that don't exist.
        """
        for key in keys:
            key_parts = key.split('_')
            if len(key_parts) > 1:
                try:
                    if int(key_parts[1]) > self.n_sources:
                        raise KeyError(
                            '{0} is defined but there are only '.format(key) +
                            '{0} sources.'.format(self.n_sources))
                except ValueError:
                    pass

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
        # If parallax is defined, then both components must be set:
        if ('pi_E_N' in keys) != ('pi_E_E' in keys):
            raise KeyError(
                'You have to define either both or none of (pi_E_N, pi_E_E).')

        # t_0_par makes sense only when parallax is defined.
        if 't_0_par' in keys and not self._type['parallax']:
            raise KeyError('t_0_par makes sense only when parallax is defined')

        # Parallax needs reference epoch:
        if 'pi_E_N' in keys:
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

        # If s_z or ds_z_dt is defined, ds_dt and dalpha_dt must be defined
        if ('s_z' in keys) and ('ds_z_dt' in keys):
            raise KeyError('Full Keplerian motion requires either s_z or ' +
                           'ds_z_dt, not both.')
        if ('s_z' in keys) or ('ds_z_dt' in keys):
            if ('ds_dt' not in keys) or ('dalpha_dt' not in keys):
                raise KeyError(
                    'Full Keplerian motion (s_z or ds_z_dt) requires ' +
                    'both ds_dt and dalpha_dt.' +
                    '\nNote that you can set either of them to 0.')
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
        required = ('period semimajor_axis inclination Omega_node argument_of_latitude_reference').split()
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

        Also, check that all values are scalars.
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
            if not np.isscalar(value) or isinstance(value, str):
                msg = "{:} must be a scalar: {:}, {:}"
                raise TypeError(msg.format(key, value, type(value)))

        for name in ['x_caustic_in', 'x_caustic_out']:
            if name in parameters.keys():
                if parameters[name] < 0. or parameters[name] > 1.:
                    msg = "Parameter {:} has to be in [0, 1] range, not {:}"
                    raise ValueError(msg.format(name, parameters[name]))

        for name in ['q']:
            if name in parameters.keys():
                if parameters[name] <= 0.:
                    msg = "Parameter {:} has to be larger than 0, not {:}"
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

    def _update_sources(self, parameter, value):
        """
        For multi-source models, update the values for all sources.
        Note that pi_E_N and pi_E_E are changed separately.
        """
        if self.n_sources == 1:
            if self.is_xallarap:
                self._update_sources_xallarap_reference()

            return

        for i in range(self.n_sources):
            source = self.__getattr__('_source_{0}_parameters'.format(i+1))
            try:
                source.__getattr__(parameter)
                source.__setattr__(parameter, value)
            except KeyError:
                continue

        if self.is_xallarap:
            if parameter == 'q_source':
                value_ = self.parameters['xi_semimajor_axis'] / value
                setattr(self._source_2_parameters, 'xi_semimajor_axis', value_)
            elif parameter == 'xi_semimajor_axis':
                value /= self.parameters['q_source']
                setattr(self._source_2_parameters, parameter, value)
            elif parameter == 'xi_argument_of_latitude_reference':
                value += 180.
                setattr(self._source_2_parameters, parameter, value)

            self._update_sources_xallarap_reference()

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

    def _get_s_z_or_ds_z_dt(self, s_z=None, ds_z_dt=None):
        """
        Calculates s_z or ds_z_dt when the other is given.
        """
        conv_factor = -self.parameters['ds_dt'] * self.parameters['s']

        if s_z is not None:
            return conv_factor / s_z
        elif ds_z_dt is not None:
            return conv_factor / ds_z_dt
        else:
            raise KeyError('s_z and ds_z_dt cannot be both None.')

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
                u_0_quantity = self.parameters['t_eff'] / self.parameters['t_E']
                return u_0_quantity
            except KeyError:
                raise AttributeError('u_0 is not defined for these parameters: {0}'.format(self.parameters.keys()))

    @u_0.setter
    def u_0(self, new_u_0):
        if self._type['Cassan08']:
            raise ValueError('u_0 cannot be set for model using ' +
                             'Cassan (2008) parameterization')
        if 'u_0' in self.parameters.keys():
            self.parameters['u_0'] = new_u_0
            self._update_sources('u_0', new_u_0)
        else:
            raise AttributeError('u_0 is not a parameter of this model.')

    @property
    def t_star(self):
        """
        *float*

        t_star = rho * t_E = source radius crossing time in days
        """
        if 't_star' in self.parameters.keys():
            return self.parameters['t_star']
        elif ('rho' in self.parameters.keys() and self._type['Cassan08']):
            return self.rho * self.t_E
        else:
            try:
                return (self.parameters['t_E'] * self.parameters['rho'])
            except KeyError:
                raise AttributeError(
                    't_star is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star.setter
    def t_star(self, new_t_star):
        if 't_star' in self.parameters.keys():
            self.parameters['t_star'] = new_t_star
            self._update_sources('t_star', new_t_star)
        else:
            raise AttributeError('t_star is not a parameter of this model.')

        if new_t_star < 0.:
            raise ValueError(
                'Source crossing time cannot be negative:', new_t_star)

    @property
    def t_eff(self):
        """
        *float*

        t_eff = u_0 * t_E = effective timescale in days
        """
        if 't_eff' in self.parameters.keys():
            return self.parameters['t_eff']
        else:
            try:
                return (self.parameters['t_E'] * self.parameters['u_0'])
            except KeyError:
                raise AttributeError(
                    't_eff is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_eff.setter
    def t_eff(self, new_t_eff):
        if 't_eff' in self.parameters.keys():
            self.parameters['t_eff'] = new_t_eff
            self._update_sources('t_eff', new_t_eff)
        else:
            raise AttributeError('t_eff is not a parameter of this model.')

    @property
    def t_E(self):
        """
        *float*

        The Einstein timescale in days.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['t_E']
        if 't_E' in self.parameters.keys():
            return self.parameters['t_E']
        elif ('t_star' in self.parameters.keys() and
              'rho' in self.parameters.keys()):
            return self.t_star / self.rho
        elif ('t_eff' in self.parameters.keys() and
              'u_0' in self.parameters.keys()):
            return self.t_eff / abs(self.u_0)
        else:
            raise AttributeError("You're trying to access t_E that was not set")

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
            self.parameters['t_E'] = new_t_E
            self._update_sources('t_E', new_t_E)
        else:
            raise AttributeError('t_E is not a parameter of this model.')

    @property
    def rho(self):
        """
        *float*

        source size as a fraction of the Einstein radius
        """
        if 'rho' in self.parameters.keys():
            return self.parameters['rho']
        elif 't_star' in self.parameters.keys() and 't_E' in self.parameters.keys():
            return self.t_star / self.t_E
        elif 't_star' in self.parameters.keys() and self._type['Cassan08']:
            return self.t_star / self.t_E
        else:
            raise AttributeError("rho is not defined and cannot be calculated")

    @rho.setter
    def rho(self, new_rho):
        if 'rho' in self.parameters.keys():
            if new_rho < 0.:
                raise ValueError('source size (rho) cannot be negative')
            self.parameters['rho'] = new_rho
            self._update_sources('rho', new_rho)
        else:
            raise AttributeError('rho is not a parameter of this model.')

    @property
    def alpha(self):
        """
        *float*

        The angle of the source trajectory relative to the binary lens
        axis (or primary-secondary axis). Measured counterclockwise,
        i.e., according to convention advocated by
        `Skowron et al. 2011 (ApJ, 738, 87) <https://ui.adsabs.harvard.edu/abs/2011ApJ...738...87S/abstract>`_,
        but shifted by 180 deg.
        """
        if self._type['Cassan08']:
            self._get_standard_parameters_from_Cassan08()
            return self._standard_parameters['alpha']

        return self.parameters['alpha']

    @alpha.setter
    def alpha(self, new_alpha):
        if self._type['Cassan08']:
            raise ValueError('alpha cannot be set for model using Cassan (2008) parameterization')

        self.parameters['alpha'] = new_alpha
        self._update_sources('alpha', new_alpha)

    @property
    def q(self):
        """
        *float*

        mass ratio of the two lens components. Only 2 bodies allowed.
        """
        return self.parameters['q']

    @q.setter
    def q(self, new_q):
        if new_q <= 0.:
            raise ValueError('mass ratio q has to be larger than 0')
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
            raise AttributeError('convergence_K is not a parameter of this model.')

    @convergence_K.setter
    def convergence_K(self, new_K):
        if 'convergence_K' in self.parameters.keys():
            self.parameters['convergence_K'] = new_K
            self._update_sources('convergence_K', new_K)
        else:
            raise AttributeError('convergence_K is not a parameter of this model.')

    @property
    def shear_G(self):
        """
        *complex*

        Shear of external mass sheet.
        """
        if 'shear_G' in self.parameters.keys():
            return self.parameters['shear_G']
        else:
            raise AttributeError('shear_G is not a parameter of this model.')

    @shear_G.setter
    def shear_G(self, new_G):
        if 'shear_G' in self.parameters.keys():
            self.parameters['shear_G'] = new_G
            self._update_sources('shear_G', new_G)
        else:
            raise AttributeError('shear_G is not a parameter of this model.')

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
            raise AttributeError('pi_E_N not defined for this model')

    @pi_E_N.setter
    def pi_E_N(self, new_value):
        if 'pi_E_N' in self.parameters.keys():
            self.parameters['pi_E_N'] = new_value
            self._update_sources('pi_E_N', new_value)
        else:
            raise AttributeError('pi_E_N is not a parameter of this model.')

    @property
    def pi_E_E(self):
        """
        *float*

        The East component of the microlensing parallax vector.
        """
        if 'pi_E_E' in self.parameters.keys():
            return self.parameters['pi_E_E']
        else:
            raise AttributeError('pi_E_N not defined for this model')

    @pi_E_E.setter
    def pi_E_E(self, new_value):
        if 'pi_E_E' in self.parameters.keys():
            self.parameters['pi_E_E'] = new_value
            self._update_sources('pi_E_E', new_value)
        else:
            raise AttributeError('pi_E_E is not a parameter of this model.')

    @property
    def pi_E(self):
        """
        Not defined.

        It was used in previous versions. Use :py:attr:`~pi_E_N` and :py:attr:`~pi_E_E` instead.
        """
        raise AttributeError('pi_E is not defined. Use pi_E_N and pi_E_E instead')

    @pi_E.setter
    def pi_E(self, new_value):
        raise AttributeError('pi_E is not defined. Use pi_E_N and pi_E_E instead')

    @property
    def t_0_par(self):
        """
        *float*

        The reference time for the calculation of parallax. If not set
        explicitly, then it is assumed t_0_par = t_0. If there are multiple sources,
        t_0_1 is used.

        Note that this is a reference value and not the fitting parameter.
        It is best to fix it at the begin of calculations.
        """
        if 't_0_par' not in self.parameters.keys():
            if 't_0_kep' in self.parameters.keys():
                return self.parameters['t_0_kep']
            elif 't_0' in self.parameters.keys():
                return self.parameters['t_0']
            elif self.n_sources > 1:
                return self.t_0_1
            else:
                raise AttributeError('No valid value for setting t_0_par', self.parameters)

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
        if 'pi_E_N' in self.parameters.keys() and 'pi_E_E' in self.parameters.keys():
            pi_E_N = self.parameters['pi_E_N']
            pi_E_E = self.parameters['pi_E_E']
        else:
            raise AttributeError('pi_E not defined for this model')
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
        *float*

        Change rate of separation :py:attr:`~s` in 1/year.
        """
        return self.parameters['ds_dt']

    @ds_dt.setter
    def ds_dt(self, new_ds_dt):
        self.parameters['ds_dt'] = new_ds_dt
        self._update_sources('ds_dt', new_ds_dt)

    @property
    def dalpha_dt(self):
        """
        *float*

        Change rate of angle :py:attr:`~alpha` in deg/year.
        """
        return self.parameters['dalpha_dt']

    @dalpha_dt.setter
    def dalpha_dt(self, new_dalpha_dt):
        self.parameters['dalpha_dt'] = new_dalpha_dt
        self._update_sources('dalpha_dt', new_dalpha_dt)

    @property
    def s_z(self):
        """
        *float*

        instantaneous position in the direction perpendicular to the plane
        of the sky at time t_0_kep.
        """
        if 's_z' in self.parameters:
            return self.parameters['s_z']
        else:
            return self._get_s_z_or_ds_z_dt(ds_z_dt=self.parameters['ds_z_dt'])

    @s_z.setter
    def s_z(self, new_s_z):
        self.parameters['s_z'] = new_s_z
        self._update_sources('s_z', new_s_z)

    @property
    def ds_z_dt(self):
        """
        *float*

        Change rate of separation :py:attr:`~s_z` in 1/year.
        """
        if 'ds_z_dt' in self.parameters:
            value = self.parameters['ds_z_dt']
        else:
            value = self._get_s_z_or_ds_z_dt(s_z=self.parameters['s_z'])

        return value

    @ds_z_dt.setter
    def ds_z_dt(self, new_ds_z_dt):
        self.parameters['ds_z_dt'] = new_ds_z_dt
        self._update_sources('ds_z_dt', new_ds_z_dt)

    @property
    def t_0_kep(self):
        """
        *float*

        The reference time for the calculation of lens orbital motion.
        If not set explicitly, then it is assumed t_0_kep = t_0 (or t_0_1 for multi-source models).

        Note that this is a reference value and not the fitting parameter.
        It is best to fix it at the begin of calculations.
        """
        if 't_0_kep' not in self.parameters.keys():
            if 't_0_par' in self.parameters.keys():
                return self.parameters['t_0_par']
            elif 't_0' in self.parameters.keys():
                return self.parameters['t_0']
            elif self.n_sources > 1:
                return self.t_0_1
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
        self._update_sources('xi_period', new_value)

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
        self._update_sources('xi_semimajor_axis', new_value)

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
        self.parameters['xi_Omega_node'] = new_value
        self._update_sources('xi_Omega_node', new_value)

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
        self.parameters['xi_inclination'] = new_value
        self._update_sources('xi_inclination', new_value)

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
        self.parameters['xi_argument_of_latitude_reference'] = new_value
        self._update_sources('xi_argument_of_latitude_reference', new_value)

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
        self._update_sources('xi_eccentricity', new_value)

    @property
    def xi_omega_periapsis(self):
        """
        *float*

        The argument of periapsis of the xallarap orbit, i.e., the angle
        between the ascending node and periapsis measured in
        the direction of motion.
        The units are degrees.
        """
        return self.parameters['xi_omega_periapsis']

    @xi_omega_periapsis.setter
    def xi_omega_periapsis(self, new_value):
        self.parameters['xi_omega_periapsis'] = new_value
        self._update_sources('xi_omega_periapsis', new_value)

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
        self._update_sources('t_0_xi', new_value)

    @property
    def q_source(self):
        """
        *float*

        The mass ratio of the second and the first source.
        This is value must be positive and can be > 1.
        Defined only for xallarap binary-source models because it does not
        affect the magnification for binary-source models without xallarap.
        """
        return self.parameters['q_source']

    @q_source.setter
    def q_source(self, new_value):
        if new_value < 0.:
            raise ValueError('q_source cannot be negative')
        self.parameters['q_source'] = new_value
        self._update_sources('q_source', new_value)

    @property
    def xallarap_reference_position(self):
        """
        *np.ndarray* of shape (2, 1)

        The position of the first source at :py:attr:`~t_0_xi` relative to
        the source center of mass. It is a 2D vector that is subtracted from
        the source position along the orbit in order to calculate the shift
        caused by xallarap.
        """
        return self._xallarap_reference_position

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
            raise AttributeError('u_0_1 is not a parameter of this model.')

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
            raise AttributeError('u_0_2 is not a parameter of this model.')

    @property
    def t_star_1(self):
        """
        *float*

        t_star_1 = rho_1 * t_E_1 = source no. 1 radius crossing time in days
        """
        if 't_star_1' in self.parameters.keys():
            return self.parameters['t_star_1']
        else:
            try:
                t_E = self._source_1_parameters.parameters['t_E']
                rho = self._source_1_parameters.parameters['rho']
                return t_E * rho
            except KeyError:
                raise AttributeError(
                    't_star_1 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star_1.setter
    def t_star_1(self, new_t_star_1):
        if 't_star_1' in self.parameters.keys():
            self.parameters['t_star_1'] = new_t_star_1
            self._source_1_parameters.t_star = new_t_star_1
        else:
            raise AttributeError('t_star_1 is not a parameter of this model.')

        if new_t_star_1 < 0.:
            raise ValueError(
                'Source crossing time cannot be negative:', new_t_star_1)

    @property
    def t_star_2(self):
        """
        *float*

        t_star_2 = rho_2 * t_E_2 = source no. 2 radius crossing time in days.
        """
        if 't_star_2' in self.parameters.keys():
            return self.parameters['t_star_2']
        else:
            try:
                t_E = self._source_2_parameters.parameters['t_E']
                rho = self._source_2_parameters.parameters['rho']
                return t_E * rho
            except KeyError:
                raise AttributeError(
                    't_star_2 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @t_star_2.setter
    def t_star_2(self, new_t_star_2):
        if 't_star_2' in self.parameters.keys():
            self.parameters['t_star_2'] = new_t_star_2
            self._source_2_parameters.t_star = new_t_star_2
        else:
            raise AttributeError('t_star_2 is not a parameter of this model.')

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
            fmt = 'rho_1 is not defined for these parameters: {:}'
            raise AttributeError(fmt.format(self.parameters.keys()))

    @rho_1.setter
    def rho_1(self, new_rho_1):
        if 'rho_1' in self.parameters.keys():
            if new_rho_1 < 0.:
                raise ValueError('source size (rho_1) cannot be negative')
            self.parameters['rho_1'] = new_rho_1
            self._source_1_parameters.rho = new_rho_1
        else:
            raise AttributeError('rho_1 is not a parameter of this model.')

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
            raise AttributeError(
                'rho_2 is not defined for these parameters: {0}'.format(
                    self.parameters.keys()))

    @rho_2.setter
    def rho_2(self, new_rho_2):
        if 'rho_2' in self.parameters.keys():
            if new_rho_2 < 0.:
                raise ValueError('source size (rho_2) cannot be negative')
            self.parameters['rho_2'] = new_rho_2
            self._source_2_parameters.rho = new_rho_2
        else:
            raise AttributeError('rho_2 is not a parameter of this model.')

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

        if self._type['full keplerian motion']:
            self._set_lens_keplerian_orbit()
            sky_positions = self._lens_orbit.get_reference_plane_position(epoch)
            s_of_t = np.sqrt(np.sum(sky_positions**2, axis=0))
        else:
            s_of_t = self.s + self.ds_dt * (epoch - self.t_0_kep) / 365.25

        return s_of_t

    def get_alpha(self, epoch):
        """
        Returns the value of angle :py:attr:`~alpha` at a given epoch or
        epochs (if orbital motion parameters are set).

        Arguments :
            epoch: *float*, *list*, *np.ndarray*
                The time(s) at which to calculate :py:attr:`~alpha`.

        Returns :
            angle: *float*
                Value(s) of angle for given epochs in degrees

        """
        if 'dalpha_dt' not in self.parameters.keys():
            return self.alpha

        if isinstance(epoch, list):
            epoch = np.array(epoch)

        if self._type['full keplerian motion']:
            self._set_lens_keplerian_orbit()
            sky_positions = self._lens_orbit.get_reference_plane_position(epoch)
            alpha_of_t = self.alpha + np.arctan2(sky_positions[1, :], sky_positions[0, :]) * 180 / np.pi
        else:
            alpha_of_t = self.alpha + self.dalpha_dt * (epoch - self.t_0_kep) / 365.25

        return alpha_of_t

    @property
    def gamma_parallel(self):
        """
        *float*

        Parallel component of instantaneous velocity of the secondary
        relative to the primary in 1/year.
        It is parallel to the primary-secondary axis.
        Equals :py:attr:`~ds_dt`/:py:attr:`~s`. Cannot be set.
        """
        return self.ds_dt / self.s

    @property
    def gamma_perp(self):
        """
        *float*

        Perpendicular component of instantaneous velocity of the secondary
        relative to the primary. It is perpendicular to the primary-secondary
        axis. It has sign opposite to :py:attr:`~dalpha_dt`
        and is in rad/yr, not deg/yr. Cannot be set.
        """
        return -self.dalpha_dt * (np.pi / 180.)

    @property
    def gamma_z(self):
        """
        *float*

        Perpendicular component of instantaneous velocity of the secondary
        relative to the primary in 1/yr. It is perpendicular to the plane of the sky
        at time :py:attr:`~t_0_kep`. Equals :py:attr:`~ds_z_dt`/:py:attr:`~s`. Cannot be set.
        """
        if not self.is_keplerian():
            return None

        return self.ds_z_dt / self.s

    @property
    def gamma(self):
        """
        *float*

        Instantaneous velocity of the secondary relative to the primary in 1/year. Cannot be set.
        """
        if not self.is_keplerian():
            return (self.gamma_parallel**2 + self.gamma_perp**2)**0.5

        return (self.gamma_parallel**2 + self.gamma_perp**2 + self.gamma_z**2)**0.5

    def _set_lens_keplerian_orbit(self):
        """
        Set parameters of the lens keplerian orbit i.e. self._lens_keplerian.
        """
        position = np.array([self.s, 0, self.s_z])
        gamma = np.array([self.gamma_parallel, self.gamma_perp, self.gamma_z])
        new_input = [*list(position), *list(gamma)]
        if new_input == self._lens_keplerian_last_input:
            return

        self._lens_keplerian_last_input = new_input

        velocity = self.s * gamma  # This is in units of R_E = D_L * theta_E.
        a = np.sqrt(np.sum(position**2))
        self._lens_keplerian['semimajor_axis'] = a
        self._lens_keplerian['period'] = 2 * np.pi * a / np.sqrt(np.sum(velocity**2)) * 365.25
        h = np.cross(position, velocity)
        j = np.array([0, 1, 0])
        n = np.cross(j, h)
        self._lens_keplerian['inclination'] = np.arctan2(np.sqrt(h[0]**2+h[1]**2), h[2]) * 180. / np.pi
        self._lens_keplerian['Omega_node'] = np.arctan2(h[0], -h[1]) * 180. / np.pi
        gamma_012 = np.sqrt(np.sum(gamma**2))
        gamma_02 = np.sqrt(gamma[0]**2+gamma[2]**2)
        phi_0 = np.arctan2(-gamma[0]*gamma_012, gamma[2]*gamma_02)
        self._lens_keplerian['argument_of_latitude_reference'] = phi_0 * 180. / np.pi
#        Utils.get_angle_between_vectors(n, position)
        self._lens_keplerian['epoch_reference'] = self.t_0_kep
#        print("ORBIT:")
#        print(self._lens_keplerian)
        self._lens_orbit = Orbit(**self._lens_keplerian)

    @property
    def lens_semimajor_axis(self):
        """
        *float*

        Semi-major axis of the binary lens orbit in units of theta_E.
        """
        self._set_lens_keplerian_orbit()
        return self._lens_keplerian['semimajor_axis']

    @property
    def lens_period(self):
        """
        *float*

        Orbital period of the binary lens orbit in years.
        """
        self._set_lens_keplerian_orbit()
        return self._lens_keplerian['period'] / 365.25

    @property
    def lens_inclination(self):
        """
        *float*

        Inclination of the binary lens orbit in degrees.
        """
        self._set_lens_keplerian_orbit()
        return self._lens_keplerian['inclination']

    @property
    def lens_Omega_node(self):
        """
        *float*

        Longitude of ascending node of the binary lens orbit in degrees.
        """
        self._set_lens_keplerian_orbit()
        return self._lens_keplerian['Omega_node']

    @property
    def lens_argument_of_latitude_reference(self):
        """
        *float*

        Argument of latitude of the binary lens orbit in degrees.
        """
        self._set_lens_keplerian_orbit()
        return self._lens_keplerian['argument_of_latitude_reference']

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
                *True* if :py:attr:`~dalpha_dt` or :py:attr:`~ds_dt` are set.

        """
        return not self._type['lens 2-parameter orbital motion']

    def is_keplerian(self):
        """
        Checks if model includes keplerian orbital motion of the lenses,
        which can be either circular or elliptical.

        Returns :
            is_keplerian: *boolean*
                *True* if :py:attr:`~s_z` or :py:attr:`~ds_z_dt` are set.
        """
        return self._type['full keplerian motion']

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
