import numpy as np

from MulensModel.modelparameters import ModelParameters

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
                 t_0=0., u_0=None, t_E=0., rho=None, s=None, q=None,
                 alpha=None,
                 lens=None, source=None, mu_rel=None):
        """
        Three ways to define the model:
        parameters = a ModelParameters() object
        specify t_0, u_0, t_E (optionally: rho, s, q, alpha)
        specify physical properties: lens= a Lens() object, 
            source= a Source() object, mu_rel
        """
        self.parameters = ModelParameters(t_0=t_0, u_0=u_0, t_E=t_E)
        if parameters is not None:
            pass
        if t_0 is not None:
            self.t_0 = t_0
        if u_0 is not None:
            self.u_0 = u_0
        if t_E is not None:
            self.t_E = t_E
        if rho is not None:
            self.rho = rho
        if source is not None:
            pass
        self._magnification = None

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
    

