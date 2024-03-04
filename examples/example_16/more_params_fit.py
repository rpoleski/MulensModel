import sys
import yaml
import numpy as np

from ulens_model_fit import UlensModelFit, import_failed

class UlensModelFitVariableBaseline(UlensModelFit):
    def _set_default_parameters(self):
        """
        Extend the set of available parameters
        """
        super()._set_default_parameters()
        self._other_parameters = ['amplitude', 'baseline_period', 'baseline_t0']
        self._latex_conversion_other = {
            'amplitude': 'amp_{\\rm base}',
            'baseline_period': 'P_{\\rm base}',
            'baseline_t0': 'epoch_{\\rm base}',
            }

    def _get_ln_probability_for_other_parameters(self):
        """
        We have to define this function and it has to return a float but in
        this case the change of ln_prob is coded in _ln_like().
        """
        return 0.

    def _ln_like(self, theta):
        """
        likelihood function
        """
        self._set_model_parameters(theta)

        # changed - getting parameters:
        params = {**self._other_parameters_dict}
        if self._fixed_parameters is not None:
            params = {**params, **self._fixed_parameters}
        amp = params['amplitude']
        period = params['baseline_period']
        t0 = params['baseline_t0']

        # changed - correcting fluxes:
        orig_fluxes = []
        for dataset in self._datasets:
            orig_fluxes.append(np.copy(dataset._flux))
            dataset._flux -= amp * np.sin(2*np.pi*(dataset.time-t0)/period)
            dataset._mag = None
            dataset._err_mag = None

        chi2 = self._event.get_chi2()
        out = -0.5 * chi2

        # changed - restoring fluxes:
        for (flux, dataset) in zip(orig_fluxes, self._datasets):
            dataset._flux = flux
            dataset._mag = None
            dataset._err_mag = None

        if self._print_model:
            self._print_current_model(theta, chi2)

        if self._task == 'fit' and len(self._other_parameters_dict) > 0:
            out += self._get_ln_probability_for_other_parameters()

        return out

# https://github.com/rpoleski/MulensModel/compare/master...ex16_galactic_model

# _get_samples_for_triangle_plot
# _get_labels_for_triangle_plot
# _get_parameters
# Wrong chi2 is printed
# wrong model is plotted


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise ValueError('Exactly one argument needed - YAML file')
    if 'yaml' in import_failed:
        raise ImportError('module "yaml" could not be imported :(')

    input_file = sys.argv[1]

    with open(input_file, 'r') as data:
        settings = yaml.safe_load(data)

    ulens_model_fit = UlensModelFitVariableBaseline(**settings)

    ulens_model_fit.run_fit()
