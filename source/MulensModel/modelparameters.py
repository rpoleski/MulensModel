from astropy import units as u

from MulensModel.mulensparallaxvector import MulensParallaxVector


class ModelParameters(object):
    """
    A class for the basic microlensing model parameters (t_0, u_0,
    t_E, rho, s, q, alpha, pi_E). Can handle point lens or binary
    lens. pi_E assumes NE coordinates (Parallel, Perpendicular
    coordinates are not supported).
    """
    def __init__(self, t_0=None, u_0=None, t_E=None, rho=None, s=None,
                 q=None, alpha=None, pi_E=None, pi_E_N=None, pi_E_E=None,
                 pi_E_ref=None):
        """
        Set up parameters for a MulensModel.Model object.
        
        Attributes:
            t_0: time of closest approach between source and lens
            u_0: impact parameter between source and lens (in Einstein radii)
            t_E: Einstein crossing time
            rho: source size as a fraction of the Einstein radius
            s: separation between primary and companion(s) (in Einstein radii)
            q: mass ratio between primary and companion(s)
            alpha: angle of source trajectory relative to binary lens axis 
            (CCW???)

            parallax vector may be defined 
            EITHER as:
               pi_E: list, tuple, or numpy.ndarray of 2 values
               pi_E_ref: defines reference system for pi_E (see
                   MulensParallaxVector)
            OR:
               pi_E_N: North component of the parallax
               pi_E_E: East component of the parallax
               
        """
        #Initialize standard parameters.
        self.t_0 = t_0
        self.u_0 = u_0
        self.t_E = t_E
        self.rho = rho
        self.s = s
        self.q = q
        self.alpha = alpha

        """
        Define the parallax if appropriate. Does not check for
        collisions (e.g. if the user specifies both pi_E and (pi_E_N,
        pi_E_E).
        """
        self._pi_E = None
        if pi_E is not None and (pi_E_N is not None or pi_E_E is not None):
            msg = 'Microlensing parallax specified in 2 ways at the same time'
            raise ValueError(msg)
        if pi_E is not None:
            self._pi_E = MulensParallaxVector(pi_E=pi_E, ref=pi_E_ref)
        if pi_E_N is not None:
            if pi_E_E is not None:
                self._pi_E = MulensParallaxVector(pi_E_1=pi_E_N, 
                                                   pi_E_2=pi_E_E, 
                                                   ref="NorthEast")
            else:
                raise AttributeError('pi_E has 2 components')

    def __repr__(self):
        """A nice way to represent a ModelParameters object as a string"""        
        #Initialize Header line
        variables = '{0:>11} {1:>9} '.format(
            "t_0 (HJD')", 'u_0')
        try:
            variables = '{0} {1:>9}'.format(
                variables, 't_E ({0})'.format(self._t_E.unit))
        except AttributeError:
            variables = '{0} {1:>9}'.format(variables, 't_E')

        #t_0 value
        try:
            values = '{0:>11.5f}'.format(self.t_0)
        except AttributeError:
            values = '{0:>11}'.format(None)
        #u_0 value
        try:
            values = '{0} {1:>9.6f}'.format(values, self._u_0)
        except AttributeError:
            values = '{0} {1:>9}'.format(values, None)

        #t_E value
        try:
            values = '{0} {1:>10.4f}'.format(values, self.t_E)
        except AttributeError:
            values = '{0} {1:>10}'.format(values,None)

        #rho value and header column
        try:
            values = '{0} {1:>7.5f}'.format(values, self._rho)
        except (AttributeError, TypeError):
            pass
        else:
            variables = '{0} {1:>7}'.format(variables, 'rho') 

        #s, q values and header columns
        try:
            variables = '{0} {1:>9} {2:>12} {3:>11}'.format(
                variables, 's', 'q', 'alpha ({0})'.format(self._alpha.unit))
            values = '{0} {1:>9.5f} {2:>12.8f} {3:>11.5f}'.format(
                values, self._s, self._q, self._alpha.value)
        except AttributeError:
            pass

        return 'Model Parameters:\n{0}\n{1}\n'.format(variables, values)

    @property
    def n_lenses(self):
        """number of objects in the lens system"""
        if self._s is None and self._q is None and self._alpha is None:
            return 1
        else:
            return 2

    @property
    def t_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self._t_0

    @t_0.setter
    def t_0(self, new_t_0):
        self._t_0 = new_t_0

    @property
    def u_0(self):
        """
        The minimum projected separation between the source
        and the lens center of mass.
        """
        if self._u_0 is not None:
            return self._u_0
        else:
            raise AttributeError('u_0 is not defined.')

    @u_0.setter
    def u_0(self, new_u_0):
        self._u_0 = new_u_0

    @property
    def t_E(self):
        """
        The Einstein timescale. An astropy.Quantity. "day" is the default unit.
        """
        return self._t_E.value 
        # Add some unit check, so that it doesn't return 1. for 1*u.year
    
    @t_E.setter
    def t_E(self, new_t_E):
        if new_t_E is None:
            self._t_E = None
            return
        if new_t_E < 0.:
            raise ValueError('Einstein timescale cannot be negative:', new_t_E)
        if isinstance(new_t_E, u.Quantity): 
            # Add a check if the unit is time unit?
            self._t_E = new_t_E
        else:
            self._t_E = new_t_E * u.day
    
    @property
    def rho(self):
        """source size as a fraction of the Einstein radius"""
        return self._rho
    
    @rho.setter
    def rho(self, new_rho):
        if new_rho is None:
            self._rho = None
            return
        if new_rho < 0.:
            raise ValueError('source size (rho) cannot be negative')
        self._rho = new_rho
    
    @property
    def alpha(self):
        """
        The angle of the source trajectory relative to the binary lens axis
        (or primary-secondary axis). Measured CW/CCW (TBD). An
        astropy.Quantity. "deg" is the default unit.
        TBD - make sure CW/CCW convention is according to Skowron+11 appendix A
        """
        return self._alpha

    @alpha.setter
    def alpha(self, new_alpha):
        if new_alpha is None:
            self._alpha = None
            return
        if isinstance(new_alpha, u.Quantity):
            self._alpha = new_alpha
        else:
            self._alpha = new_alpha * u.deg

    @property
    def q(self):
        """mass ratio of two lens components"""
        return self._q
        
    @q.setter
    def q(self, new_q):
        self._q = new_q
    
    @property
    def s(self):
        """separation of two lens components relative to Einstein ring size"""
        return self._s

    @s.setter
    def s(self, new_s):
        self._s = new_s

    @property
    def pi_E(self):
        """
        The microlens parallax vector. May be specified either
        relative to the sky ("NorthEast") or relative to the binary lens 
        axis ("ParPerp"). "NorthEast" is default.
        """
        return self._pi_E

    @pi_E.setter
    def pi_E(self, new_pi_E):
        if isinstance(new_pi_E, MulensParallaxVector):
            self._pi_E = new_pi_E
        else:
            self._pi_E = MulensParallaxVector(pi_E=new_pi_E, ref=None)

    @property
    def pi_E_N(self):
        """
        The North component of the microlens parallax vector.
        """
        return self._pi_E.vector[0]

    @pi_E_N.setter
    def pi_E_N(self, new_value):
        try:
            self._pi_E.vector[0] = new_value
        except AttributeError:
            self._pi_E = MulensParallaxVector(pi_E_1=new_value, pi_E_2=0., 
                                               ref="NorthEast")
    @property
    def pi_E_E(self):
        """
        The East component of the microlens parallax vector.
        """
        return self._pi_E.vector[1]

    @pi_E_E.setter
    def pi_E_E(self, new_value):
        try:
            self._pi_E.vector[1] = new_value
        except AttributeError:
            self._pi_E = MulensParallaxVector(pi_E_1=0., pi_E_2=new_value, 
                                               ref="NorthEast")
