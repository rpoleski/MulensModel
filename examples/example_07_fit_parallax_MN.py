"""
Fits parallax PSPL model using MultiNest technique. It finds two
separate modes in automated way.

"""
import os
import sys
import numpy as np
try:
    from pymultinest.solve import solve
    from pymultinest.analyse import Analyzer
except ImportError as err:
    print(err)
    print("\nPyMultiNest could not be imported.")
    print("Get it from: https://github.com/JohannesBuchner/PyMultiNest")
    print("and re-run the script")
    sys.exit(1)

import MulensModel as mm


class Minimizer(object):
    """
    A class used to store settings and functions that are used and
    called by MultiNest.

    Parameters :
        event: *MulensModel.event*
            Event for which chi^2 will be calculated.

        parameters_to_fit: *list*
            Names of parameters to be fitted, e.g.,  ["t_0", "u_0", "t_E"]

    """
    def __init__(self, event, parameters_to_fit):
        self.event = event
        self.parameters_to_fit = parameters_to_fit
        self._chi2_0 = None

    def set_chi2_0(self, chi2_0=None):
        """set reference value of chi2 (can help with numerical stability"""
        if chi2_0 is None:
            chi2_0 = np.sum([d.n_epochs for d in self.event.datasets])
        self._chi2_0 = chi2_0

    def set_cube(self, min_values, max_values):
        """
        remembers minimum and maximum values of model parameters so that later
        transform_cube() can be called
        """
        self._zero_points = min_values
        self._differences = max_values - min_values

    def transform_cube(self, cube):
        """transforms n-dimensional cube into space of physical quantities"""
        return self._zero_points + self._differences * cube

    def chi2(self, theta):
        """return chi^2 for a parameters theta"""
        for (i, param) in enumerate(self.parameters_to_fit):
            setattr(self.event.model.parameters, param, theta[i])
        chi2 = self.event.get_chi2()
        return chi2

    def ln_likelihood(self, theta):
        """logarithm of likelihood"""
        return -0.5 * (self.chi2(theta) - self._chi2_0)


# Read the data
file_name = os.path.join(
    mm.DATA_PATH, "photometry_files", "OB05086", "starBLG234.6.I.218982.dat")
my_data = mm.MulensData(file_name=file_name, add_2450000=True)

# Starting parameters:
coords = "18:04:45.71 -26:59:15.2"
t_0_par = 2453628.
params = {'t_0': 2453628.3, 't_0_par': t_0_par, 'u_0': 0.37, 't_E': 100.,
          'pi_E_N': 0., 'pi_E_E': 0.}
my_model = mm.Model(params, coords=coords)
my_event = mm.Event(datasets=my_data, model=my_model)

# Which parameters we want to fit?
parameters_to_fit = ["t_0", "u_0", "t_E", "pi_E_N", "pi_E_E"]
min_values = np.array([t_0_par - 5., -0.6, 80., -0.5, -0.5])
max_values = np.array([t_0_par + 5., 0.6, 120., 0.5, 0.5])

minimizer = Minimizer(my_event, parameters_to_fit)
minimizer.set_cube(min_values, max_values)
minimizer.set_chi2_0()

# Make sure the directory with detailed results exists:
dir_out = "chains/"
if not os.path.exists(dir_out):
    os.mkdir(dir_out)
file_prefix = __file__.split(".py")[0]

# Run MultiNest:
run_kwargs = {
    'LogLikelihood': minimizer.ln_likelihood,
    'Prior': minimizer.transform_cube,
    'n_dims': len(parameters_to_fit),
    'outputfiles_basename': os.path.join(dir_out, file_prefix+"_"),
    'resume': False,
    'importance_nested_sampling': False,
    'multimodal': True}
result = solve(**run_kwargs)

# Analyze results - we print each mode separately and give log-evidence:
analyzer = Analyzer(n_params=run_kwargs['n_dims'],
                    outputfiles_basename=run_kwargs['outputfiles_basename'])
modes = analyzer.get_mode_stats()['modes']
for mode in modes:
    print("\nMode {:} parameters:".format(mode['index']))
    for (p, m, s) in zip(parameters_to_fit, mode['mean'], mode['sigma']):
        print("{:6} = {:.4f} +- {:.4f}".format(p, m, s))
    print("local log-evidence: {:.3f} +- {:.3f}".format(
            mode["local log-evidence"], mode["local log-evidence error"]))
