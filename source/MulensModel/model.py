
class Model(object):
    def __init__(self,parameters=None, lens=None,source=None):
        pass

    @property
    def t_0(self):
        return self.parameters._t_0

    @property
    def u_0(self):
        return self.parameters._u_0

    @property
    def t_E(self):
        return self.parameters._t_E

class _Parameters(object):
    def __init__(self,t_0=0., u_0=0., t_E=0., rho=None):

