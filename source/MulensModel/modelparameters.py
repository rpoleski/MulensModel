from astropy import units as u
import numpy as np


# For definition of class ModelParameters see below.

# Different parameter sets. Any parameters that may be given as
# 'basic' should be a list. Parameters that may be 'optional' should
# be a list of length 2. The second item will only be printed if the
# effect is included in the 'optional' list (see _get_effect_strings()).
_valid_parameters = {
    'point lens': ['t_0, u_0, t_E'],
    'point lens alt': 'alternate: t_eff may be substituted for u_0 or t_E',
    'binary lens': ['s, q, alpha'],
    'finite source': ['rho', '(for finite source effects)'],
    'finite source alt': 'alternate: t_star may be substituted for t_E or rho',
    'parallax': ['pi_E OR pi_E_N, pi_E_E', '(for parallax)'],
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
    q, alpha], alternate = [teff, tstar], and optional = [pi_E or
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
    if len(args) == 0:
        _print_all()
    else:
        components = _get_effect_strings(*args)
        header = '---------\n{0} parameters:'.format(args[0])
        _print_parameters(header, components)


class ModelParameters(object):
    """
    A class for the basic microlensing model parameters (t_0, u_0,
    t_E, rho, s, q, alpha, pi_E). Can handle point lens or binary
    lens. The pi_E assumes NE coordinates (Parallel, Perpendicular
    coordinates are not supported).

    Arguments :
        parameters: *dictionary*

            A dictionary of parameters and their values. See
            :py:func:`which_parameters()` for valid parameter combinations.

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

        self._check_valid_combination(parameters.keys())
        self._set_parameters(parameters)

    def __repr__(self):
        """A nice way to represent a ModelParameters object as a string"""
        variables, values = '', ''
        if 't_0' in self.parameters.keys():
            variables += '{0:>13} '.format('t_0 (HJD)')
            values += '{0:>13.5f} '.format(self.t_0)

        if 'u_0' in self.parameters.keys():
            variables += '{0:>9} '.format('u_0')
            values += '{0:>9.6f} '.format(self.u_0)

        if 't_eff' in self.parameters.keys():
            variables += '{0:>10} '.format('t_eff (d)')
            values += '{0:>10.6f} '.format(self.t_eff)

        if 't_E' in self.parameters.keys():
            variables += '{0:>10} '.format('t_E (d)')
            values += '{0:>10.4f} '.format(self.t_E)

        if 'rho' in self.parameters.keys():
            variables += '{0:>7} '.format('rho')
            values += '{0:>7.5f} '.format(self.rho)

        if 't_star' in self.parameters.keys():
            variables += '{0:>10} '.format('t_star (d)')
            values += '{0:>10.6f} '.format(self.t_star)

        if ('pi_E' in self.parameters.keys() or
                'pi_E_N' in self.parameters.keys()):
            variables += '{0:>9} {1:>9} '.format('pi_E_N', 'pi_E_E')
            values += '{0:>9.5f} {1:>9.5f} '.format(self.pi_E_N, self.pi_E_E)

        if ('s' in self.parameters.keys() or
                'q' in self.parameters.keys() or
                'alpha' in self.parameters.keys()):
            variables += '{0:>9} {1:>12} {2:>11} '.format(
                's', 'q', 'alpha ({0})'.format(self.alpha.unit))
            values += '{0:>9.5f} {1:>12.8f} {2:>11.5f} '.format(
                self.s, self.q, self.alpha.value)
            if 'ds_dt' in self.parameters.keys():
                variables += '{0:>11} {1:>18} '.format(
                    'ds/dt (/yr)', 'dalpha/dt (deg/yr)')
                values += '{0:11.5f} {1:>18.5f} '.format(
                    self.ds_dt, self.dalpha_dt)

        return '{0}\n{1}\n'.format(variables, values)

    def _check_valid_combination(self, keys):
        """
        Check that the user hasn't over-defined the ModelParameters.

        """
        # ***Check minimal parameters for a model are defined (Not
        # Implemented)***

        # If s, q, and alpha must all be defined if one is defined
        if ('s' in keys) or ('q' in keys) or ('alpha' in keys):
            if ('s' not in keys) or ('q' not in keys) or ('alpha' not in keys):
                raise KeyError(
                    'A binary model requires all three of (s, q, alpha).')

        # Cannot define all 3 parameters for 2 observables
        if ('t_E' in keys) and ('rho' in keys) and ('t_star' in keys):
            raise KeyError('Only 1 or 2 of (t_E, rho, t_star) may be defined.')

        if ('t_E' in keys) and ('u_0' in keys) and ('t_eff' in keys):
            raise KeyError('Only 1 or 2 of (u_0, t_E, t_eff) may be defined.')

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
        if 't_0_par' in keys:
            if 'pi_E' not in keys and 'pi_E_N' not in keys:
                raise KeyError(
                    't_0_par makes sense only when parallax is defined.')

        # Parallax needs reference epoch:
        if 'pi_E' in keys or 'pi_E_N' in keys:
            if 't_0' not in keys and 't_0_par' not in keys:
                raise KeyError(
                    'Parallax is defined, hence either t_0 or t_0_par has ' +
                    'to be set.')

        # If ds_dt is defined, dalpha_dt must be defined
        if ('ds_dt' in keys) or ('dalpha_dt' in keys):
            if ('ds_dt' not in keys) or ('dalpha_dt' not in keys):
                raise KeyError(
                    'Lens orbital motion requires both ds_dt and dalpha_dt.' +
                    '\nNote that you can set either of them to 0.')
        # If orbital motion is defined, then we need binary lens.
            if (('s' not in keys) or ('q' not in keys) or
                  ('alpha' not in keys)):
                raise KeyError(
                    'Lens orbital motion requires >2 bodies (s, q, alpha).')
        # If orbital motion is defined, then reference epoch has to be set.
            if 't_0' not in keys and 't_0_kep' not in keys:
                raise KeyError('Orbital motion requires reference epoch, ' +
                    'i.e., t_0 or t_0_kep')

        # t_0_kep makes sense only when orbital motion is defined.
        if 't_0_kep' in keys:
            if 'ds_dt' not in keys or 'dalpha_dt' not in keys:
                raise KeyError(
                    't_0_kep makes sense only when orbital motion is defined.')

        # Make sure user does not set the gamma parameters.
        if 'gamma' in keys or 'gamma_perp' in keys or 'gamma_parallel' in keys:
            raise KeyError('You cannot set gamma, gamma_perp, ' +
                'or gamma_parallel. These are derived parameters. ' +
                'You can set ds_dt and dalpha_dt instead.')

    def _check_valid_parameter_values(self, parameters):
        """
        Prevent user from setting negative (unphysical) values for
        t_E, t_star, rho.
        """
        names = ['t_E', 't_star', 'rho']
        full_names = {
            't_E': 'Einstein timescale',
            't_star': 'Source crossing time', 'rho': 'Source size'}

        for name in names:
            if name in parameters.keys():
                if parameters[name] < 0.:
                    raise ValueError("{:} cannot be negative: {:}".format(
                            full_names[name], parameters[name]))

    def _set_parameters(self, parameters):
        self._check_valid_parameter_values(parameters)
        self.parameters = dict(parameters)

    @property
    def t_0(self):
        """
        *float*

        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self.parameters['t_0']

    @t_0.setter
    def t_0(self, new_t_0):
        self.parameters['t_0'] = new_t_0

    @property
    def u_0(self):
        """
        *float*

        The minimum projected separation between the source
        and the lens center of mass.
        """
        if 'u_0' in self.parameters.keys():
            return self.parameters['u_0']
        else:
            try:
                return self.parameters['t_eff'] / self.parameters['t_E']
            except KeyError:
                raise AttributeError(
                    'u_0 is not defined for these parameters: {0}'.format(
                        self.parameters.keys()))

    @u_0.setter
    def u_0(self, new_u_0):
        if 'u_0' in self.parameters.keys():
            self.parameters['u_0'] = new_u_0
        else:
            raise KeyError('u_0 is not a parameter of this model.')

    @property
    def t_star(self):
        """
        *float*

        t_star = rho * tE = source radius crossing time

        "day" is the default unit. Regardless of input value, returns
        value with units of u.day. May be set as a *float* --> assumes
        units of degrees.

        Returns:
            *float* value in days.
        """
        if 't_star' in self.parameters.keys():
            self._check_time_quantity('t_star')
            return self.parameters['t_star'].to(u.day).value
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

        "day" is the default unit. Regardless of input value, returns
        value with units of u.day. May be set as a *float* --> assumes
        units of degrees.

        Returns:
            *float* value in days.
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
        else:
            raise KeyError('t_eff is not a parameter of this model.')

    @property
    def t_E(self):
        """
        *astropy.Quantity*

        The Einstein timescale. "day" is the default unit. Regardless
        of input value, returns value with units of u.day. May be set
        as a *float* --> assumes units of degrees.
        """
        if 't_E' in self.parameters.keys():
            self._check_time_quantity('t_E')
            return self.parameters['t_E'].to(u.day).value
        elif ('t_star' in self.parameters.keys() and
              'rho' in self.parameters.keys()):
            return self.t_star/self.rho
        elif ('t_eff' in self.parameters.keys() and
              'u_0' in self.parameters.keys()):
            return self.t_eff/self.u_0
        else:
            raise KeyError("You're trying to access t_E that was not set")

    @t_E.setter
    def t_E(self, new_t_E):
        if new_t_E is None:
            raise ValueError('Must provide a value')

        if new_t_E < 0.:
            raise ValueError('Einstein timescale cannot be negative:', new_t_E)

        if 't_E' in self.parameters.keys():
            self._set_time_quantity('t_E', new_t_E)
        else:
            raise KeyError('t_E is not a parameter of this model.')

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
        if not isinstance(self.parameters[key], u.Quantity):
            self._set_time_quantity(key, self.parameters[key])

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
            return self.t_star/self.t_E
        else:
            return None

    @rho.setter
    def rho(self, new_rho):
        if 'rho' in self.parameters.keys():
            if new_rho < 0.:
                raise ValueError('source size (rho) cannot be negative')
            self.parameters['rho'] = new_rho
        else:
            raise KeyError('rho is not a parameter of this model.')

    @property
    def alpha(self):
        """
        *astropy.Quantity*

        The angle of the source trajectory relative to the binary lens
        axis (or primary-secondary axis). Measured counterclockwise,
        i.e., according to convention advocated by `Skowron et
        al. 2011 (ApJ, 738, 87)
        <http://adsabs.harvard.edu/abs/2011ApJ...738...87S>`_.  May be
        set as a *float* --> assumes "deg" is the default unit.
        Regardless of input value, returns value in degrees.
        """
        if not isinstance(self.parameters['alpha'], u.Quantity):
            self.parameters['alpha'] = self.parameters['alpha'] * u.deg

        return self.parameters['alpha'].to(u.deg)

    @alpha.setter
    def alpha(self, new_alpha):
        if isinstance(new_alpha, u.Quantity):
            self.parameters['alpha'] = new_alpha
        else:
            self.parameters['alpha'] = new_alpha * u.deg

    @property
    def q(self):
        """
        *float*

        mass ratio of the two lens components. Only 2 bodies allowed.
        """
        if isinstance(self.parameters['q'], (list, np.ndarray)):
            self.parameters['q'] = self.parameters['q'][0]
        return self.parameters['q']

    @q.setter
    def q(self, new_q):
        self.parameters['q'] = new_q

    @property
    def s(self):
        """
        *float*

        separation of the two lens components relative to Einstein ring size
        """
        if isinstance(self.parameters['s'], (list, np.ndarray)):
            self.parameters['s'] = self.parameters['s'][0]
        return self.parameters['s']

    @s.setter
    def s(self, new_s):
        self.parameters['s'] = new_s

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
            else:
                raise TypeError('pi_E is a 2D vector. It must have length 2.')

        elif ('pi_E_N' in self.parameters.keys() and
              'pi_E_E' in self.parameters.keys()):
            self.parameters['pi_E_N'] = new_pi_E[0]
            self.parameters['pi_E_E'] = new_pi_E[1]
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
        elif 'pi_E' in self.parameters.keys():
            self.parameters['pi_E'][0] = new_value
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
        elif 'pi_E' in self.parameters.keys():
            self.parameters['pi_E'][1] = new_value
        else:
            raise KeyError('pi_E_E is not a parameter of this model.')

    @property
    def t_0_par(self):
        """
        *float*

        The reference time for the calculation of parallax. If not set
        explicitly, set t_0_par = t_0.
        """
        if 't_0_par' not in self.parameters.keys():
            return self.parameters['t_0']
        else:
            return self.parameters['t_0_par']

    @t_0_par.setter
    def t_0_par(self, new_t_0_par):
        self.parameters['t_0_par'] = new_t_0_par

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
    def ds_dt(self):
        """
        *astropy.Quantity*

        Change rate of separation py:attr:`~s` in 1/year. Can be set as
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

    @property
    def dalpha_dt(self):
        """
        *astropy.Quantity*

        Change rate of angle py:attr:`~alpha` in deg/year. Can be set as 
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

    @property
    def t_0_kep(self):
        """
        *float*

        The reference time for the calculation of parallax. If not set
        explicitly, assumes t_0_kep = t_0.
        """
        if 't_0_kep' not in self.parameters.keys():
            return self.parameters['t_0']
        else:
            return self.parameters['t_0_kep']

    @t_0_kep.setter
    def t_0_kep(self, new):
        self.parameters['t_0_kep'] = new

    def get_s(self, epoch):
        """
        Returns the value of separation py:attr:`~s` at a given epoch or
        epochs (if orbital motion parameters are set).

        Arguments :
            epoch: *float*, *list*, *np.ndarray*
                The time(s) at which to calculate py:attr:`~s`.

        Returns :
            separation: *float* or *np.ndarray*
                Value(s) of separation for given epochs.

        """
        if 'ds_dt' not in self.parameters.keys():
            return self.s

        if isinstance(epoch, list):
            epoch = np.array(epoch)

        s_of_t = (self.s + self.ds_dt * (epoch - self.t_0_kep) * u.d)

        return s_of_t.value

    def get_alpha(self, epoch):
        """
        Returns the value of angle py:attr:`~alpha` at a given epoch or
        epochs (if orbital motion parameters are set).

        Arguments :
            epoch: *float*, *list*, *np.ndarray*
                The time(s) at which to calculate py:attr:`~alpha`.

        Returns :
            separation: *astropy.Quantity*
                Value(s) of angle for given epochs in degrees

        """
        if 'dalpha_dt' not in self.parameters.keys():
            return self.alpha

        if isinstance(epoch, list):
            epoch = np.array(epoch)

        alpha_of_t = (self.alpha + self.dalpha_dt * (epoch - self.t_0_kep)*u.d)
        
        return alpha_of_t.to(u.deg)

    @property
    def gamma_parallel(self):
        """
        *astropy.Quantity*

        Parallel component of instantaneous velocity of the secondary
        relative to the primary in 1/year.
        It is parallel to the primary-secondary axis.
        Equals py:attr:`~ds_dt`/py:attr:`~s`. Cannot be set.
        """
        return self.ds_dt / self.s

    @property
    def gamma_perp(self):
        """
        *astropy.Quantity*

        Perpendicular component of instantaneous velocity of the secondary
        relative to the primary. It is perpendicular to the primary-secondary
        axis. It has sign opposite to py:attr:`~dalpha_dt`
        and is in rad/yr, not deg/yr. Cannot be set.
        """
        return -self.dalpha_dt.to(u.rad/u.yr)

    @property
    def gamma(self):
        """
        *astropy.Quantity*

        Instantaneous velocity of the secondary relative to the primary in
        1/year. Cannot be set.
        """
        gamma_perp = (self.gamma_perp / u.rad).to(1/u.yr)
        return (self.gamma_parallel**2 + gamma_perp**2)**0.5

    def is_static(self):
        """
        Checks if model is static, i.e., orbital motion parameters are not set.

        Returns :
            is_static: *boolean*
                *True* if *dalpha_dt* or *ds_dt* are set.
        
        """
        if ('dalpha_dt' in self.parameters.keys() or
                'ds_dt' in self.parameters.keys()):
            return False
        else:
            return True

    @property
    def n_lenses(self):
        """
        *int*

        number of objects in the lens system
        """
        if (('s' not in self.parameters.keys()) and
                ('q' not in self.parameters.keys()) and
                ('alpha' not in self.parameters.keys())):
            return 1
        else:
            return 2

    def as_dict(self):
        """
        Give parameters as a dict.

        Returns :
            dictionary: *dict*
                The dictionary of model parameters.
        """
        return self.parameters
