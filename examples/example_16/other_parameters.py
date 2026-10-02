import sys
import yaml
import numpy as np

from ulens_model_fit import UlensModelFit

class MyUlensModelFit(UlensModelFit):
    """
    Redefines UlensModelFit but adds D_L to the fitted parameters.
    """
    def _set_default_user_and_other_parameters(self):
        self._other_parameters = ['D_L']
        self._latex_conversion_other = {'D_L': 'D_{L}'}
        self._check_if_DS_in_extras()


    def _get_ln_probability_for_other_parameters(self):
        source_distance = self._add_source_distance()
        out = self._get_ln_normal(source_distance, 1, 8.5)
        return out

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = MyUlensModelFit(**settings)

    ulens_model_fit.run_fit()
