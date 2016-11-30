import numpy as np
from astropy import units as u

from MulensModel.mulenstime import MulensTime
from MulensModel.mulensparallaxvector import MulensParallaxVector

class ModelParameters(object):
    """
    A class for the basic microlensing model parameters (t_0, u_0,
    t_E, rho, s, q, alpha, pi_E). Can handle point lens or binary
    lens. pi_E assumes NE coordinates (Parallel, Perpendicular
    coordinates are not supported).
    """
    def __init__(self, t_0=0., u_0=None, t_E=1., rho=None, s=None,
                 q=None, alpha=None, pi_E=None, pi_E_N=None, pi_E_E=None,
                 pi_E_ref=None):
        self.t_0 = t_0
        if u_0 is not None:
            self._u_0 = u_0
        self.t_E = t_E
        if rho is not None:
            self._rho = rho
        if s is not None:
            self._s = s
        if q is not None:
            self._q = q
        if alpha is not None:
            self.alpha = alpha

        """
        Define the parallax if appropriate. Does not check for
        collisions (e.g. if the user specifies both pi_E and (pi_E_N,
        pi_E_E).
        """
        if pi_E is not None and (pi_E_N is not None or pi_E_E is not None):
            msg = 'Microlensing parallax specified in 2 ways at the same time'
            raise ValueError(msg)
        if pi_E is not None:
            self._pi_E = MulensParallaxVector(pi_E, ref=pi_E_ref)
        if pi_E_N is not None:
            if pi_E_E is not None:
                self._pi_E = MulensParallaxVector(pi_E_1=pi_E_N, 
                                                   pi_E_2=pi_E_E, 
                                                   ref="NorthEast")
            else:
                raise AttributeError('pi_E has 2 components')

    def __repr__(self):
        variables = '{0:>11} {1:>9} {2:>10}'.format(
            "t_0 (HJD')", 'u_0', 
            't_E ({0})'.format(self._t_E.unit))
        try:
            values = '{0:>11.5f}'.format(self.t_0)
        except AttributeError:
            values = '{0:>11}'.format(None)
        try:
            values = '{0} {1:>9.6f}'.format(values, self._u_0)
        except AttributeError:
            values = '{0} {1:>9}'.format(values, None)
        try:
            values = '{0} {1:>10.4f}'.format(values,self.t_E)
        except AttributeError:
            values = '{0} {1:>10}'.format(values,None)
        try:
            variables = '{0} {1:>7}'.format(variables, 'rho')
            values = '{0} {1:>7.5f}'.format(values, self._rho)
        except AttributeError:
            pass
        try:
            variables = '{0} {1:>9} {2:>12} {3:>11}'.format(
                variables, 's', 'q', 'alpha ({0})'.format(self._alpha.unit))
            values = '{0} {1:>9.5f} {2:>12.8f} {3:>11.5f}'.format(
                values, self._s, self._q, self._alpha.value)
        except AttributeError:
            pass
        return 'ModelParameters:\n{0}\n{1}\n'.format(variables, values)


    @property
    def t_0(self):
        """
        The time of minimum projected separation between the source
        and the lens center of mass.
        """
        return self._t_0.time

    @t_0.setter
    def t_0(self, new_t_0):
        if isinstance(new_t_0, MulensTime):
            self._t_0 = new_t_0
        else:
            self._t_0 = MulensTime(new_t_0, date_fmt="hjdprime")

    @property
    def u_0(self):
        """
        The time of minimum projected separation between the source
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
        if new_t_E < 0.:
            raise ValueError('Einstein timescale cannot be negaitve')
        if isinstance(new_t_E, u.Quantity): 
            # Add a check if the unit is time unit?
            self._t_E = new_t_E
        else:
            self._t_E = new_t_E * u.day
    
    @property
    def alpha(self):
        """
        The angle of the source trajectory relative to the binary axis
        (or primary-secondary axis). Measured CW/CCW (TBD). An
        astropy.Quantity. "deg" is the default unit.
        TBD - make sure CW/CCW convention is according to Skowron+11 appendix A
        """
        return self._alpha

    @alpha.setter
    def alpha(self, new_alpha):
        if isinstance(new_alpha, u.Quantity):
            self._alpha = new_alpha
        else:
            self._alpha = new_alpha * u.deg

    @property
    def pi_E(self):
        """
        The microlens parallax vector. May be specified either
        relative to the sky ("NorthEast") or relative to the binary
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

if __name__ == "__main__":
    print('test __repr__ for ModelParameters()')

    params_1 = ModelParameters(t_0=7620., u_0=0.001, t_E=23.*u.day)
    print(params_1)

    params_2 = ModelParameters(t_0=7600., u_0=0.3, t_E=423.*u.day,
                               rho=0.001, s=1.4, q=0.002, alpha=0.4*u.rad)
    print(params_2)
    params_3 = ModelParameters(t_0=6100., u_0=0.24, t_E=0.6,
                               rho=0.0006, s=0.1, q=0.4, alpha=12.)
    print(params_3)

