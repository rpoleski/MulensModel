import sys
import yaml
import numpy as np

from ulens_model_fit import UlensModelFit


def get_my_parameters(parameters):
    """
    Transformation of parameters for WIDE orbit planets:
    (t_0_pl, u_0_pl, t_E_pl) ==> (s, q, alpha)
    Parameters t_0, u_0, and t_E are required in input.
    """
    t_E_ratio = parameters['t_E_pl'] / parameters['t_E']
    u = parameters['u_0'] + parameters['u_0_pl'] * t_E_ratio
    tau = (parameters['t_0_pl'] - parameters['t_0']) / parameters['t_E']
    ss = np.sqrt(u**2 + tau**2)
    parameters['s'] = 0.5 * (ss + np.sqrt(ss**2+4.))
    parameters['q'] = t_E_ratio**2
    parameters['alpha'] = np.arcsin(u/ss) * 180 / np.pi
    for p in ['t_0_pl', 'u_0_pl', 't_E_pl']:
        parameters.pop(p)
    
    return parameters


class MyUlensModelFit(UlensModelFit):
    """
    Redefines UlensModelFit but replaces (s,q,alpha) parameters with (t_0_pl,u_0_pl,t_E_pl).
    """
    def _set_default_user_and_other_parameters(self):
        self._user_parameters = ['t_0_pl', 'u_0_pl', 't_E_pl']
        self._latex_conversion_user = {'t_0_pl': 't_{0, pl}', 'u_0_pl': 'u_{0, pl}', 't_E_pl': 't_{{\\rm E}, pl}'}

    def _transform_parameters(self, p):
        return get_my_parameters(p)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = MyUlensModelFit(**settings)

    ulens_model_fit.run_fit()

