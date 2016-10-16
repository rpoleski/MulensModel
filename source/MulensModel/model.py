import numpy as np
from astropy import units as u

from MulensModel.mulenstime import MulensTime

class ModelParameters(object):
    """
    A class for the basic microlensing model parameters (t_0, u_0,
    t_E, rho, s, q, alpha). Can handle point lens or binary lens. 
    """
    def __init__(self, t_0=0., u_0=0., t_E=0., rho=None, s=None,
                 q=None, alpha=None):
        self.t_0 = t_0
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

    def __repr__(self):
        variables = '{0:>11} {1:>9} {2:>10}'.format(
            "t_0 (HJD')", 'u_0', 
            't_E ({0})'.format(self._t_E.unit))
        values = '{0:>11.5f} {1:>9.6f} {2:>10.4f}'.format(
            self.t_0, self._u_0, self.t_E)
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
    def t_0(self, new_t_0, date_fmt="hjdprime"):
        if type(new_t_0) is MulensTime:
            self._t_0 = new_t_0
        else:
            self._t_0 = MulensTime(new_t_0, date_fmt=date_fmt)

    @property
    def t_E(self):
        """
        The Einstein timescale. An astropy.Quantity. "day" is the default unit.
        """
        return self._t_E.value
    
    @t_E.setter
    def t_E(self, new_t_E):
        if type(new_t_E) is u.Quantity:
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
        if type(new_alpha) is u.Quantity:
            self._alpha = new_alpha
        else:
            self._alpha = new_alpha * u.deg

class Model(object):
    """
    Caveats:
    1. Does not currently have self-consistency checks: e.g. it is
    possible to define s, a_proj, and a source distance that are not
    self-consistent. Under these circumstances, the behavior may be
    unpredictable.
    2. Does not permit parallax
    """
    def __init__(self, parameters=None,
                 t_0=None, u_0=None, t_E=None, rho=None, s=None, q=None,
                 alpha=None,
                 lens=None, source=None, mu_rel=None):
        """
        Three ways to define the model:
        parameters = a ModelParameters() object
        specify t0, u0, tE (optionally: rho, s, q, alpha)
        specify physical properties: lens= a Lens() object, 
            source= a Source() object, mu_rel
        """
        if parameters is not None:
            pass
        elif t_0 is not None:
            pass
        elif source is not None:
            pass
        else:
            raise TypeError('Not a valid model definiion')
        self._magnification = None
        self.parameters = ModelParameters()

    @property
    def t_0(self):
        return self.parameters.t_0

    @t_0.setter
    def t_0(self, value):
        self.parameters.t_0 = value

    @property
    def u_0(self):
        return self.parameters._u_0
    
    @u_0.setter
    def u_0(self, value):
        self.parameters._u_0 = value

    @property
    def t_E(self):
        return self.parameters.t_E

    @t_E.setter
    def t_E(self, value):
        self.parameters.t_E = value

    @property
    def magnification(self):
        """a list of magnifications calculated for every dataset time vector"""
        if self._magnification is not None:
            return self._magnification
        self._magnification = []
        for dataset in self._datasets:
            time_diff = (dataset.time - self.t_0) / self.t_E
            u2 = self.u_0 * self.u_0 + time_diff * time_diff
            u = np.sqrt(u2)
            self._magnification.append((u2 + 2.) / (u * np.sqrt(u2 + 4.)))
        return self._magnification

    @magnification.setter
    def magnification(self, new_value):
        self._magnification = new_value

    def set_datasets(self, datasets):
        """set _datasets property"""
        self._datasets = datasets

    @property
    def lens(self):
        pass

    @lens.setter
    def lens(self, new_lens):
        pass

    @property
    def source(self):
        pass

    @source.setter
    def source(self, new_lens):
        pass
    

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

