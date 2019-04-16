import os
import numpy as np
import scipy.optimize as op
from matplotlib import pyplot as plt
from matplotlib import gridspec

import MulensModel as MM


def chi2_fun(theta, event, parameters_to_fit):
    """
    Calculate chi^2 for *event* with *parameters_to_fit*
    set to corresponding *theta* values.
    """
    for (theta_, param) in zip(theta, parameters_to_fit):
        setattr(event.model.parameters, param, theta_)
    return event.get_chi2()


# Read in the data file.
# Maybe you have to change the path to this file.
SAMPLE_FILE_01 = os.path.join(
    MM.MODULE_PATH, "data", "photometry_files", "OB08092",
    "phot_ob08092_O4.dat")
data = MM.MulensData(file_name=SAMPLE_FILE_01)

# Initialize the fit.
parameters = ["t_0", "u_0", "t_E"]
initial_guess = [5380., 0.5, 20.]
model = MM.Model(dict(zip(parameters, initial_guess)))
event = MM.Event(datasets=data, model=model)

# Find the best-fit parameters.
result = op.minimize(
    chi2_fun, x0=initial_guess, args=(event, parameters),
    method='Nelder-Mead')

# Print results.
print("{:.5f} {:.4f} {:.3f}".format(*result.x))
(fit_t_0, fit_u_0, fit_t_E) = result.x
chi2 = chi2_fun(result.x, event, parameters)
print("{:.3f}".format(chi2))

# Plot results.
best = [event.best_chi2_parameters[p] for p in parameters]
chi2_fun(best, event, parameters) # This allows plotting of the best model.
grid = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
plt.subplot(grid[0])
event.plot_data()
event.plot_model(t_start=5340., t_stop=5420.)
plt.xlim(5340., 5420.)
plt.subplot(grid[1])
event.plot_residuals()
plt.xlim(5340., 5420.)
plt.show()
