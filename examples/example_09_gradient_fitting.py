"""
Fit point source point lens model to OB08092 using Newton-CG method from scipy.
This method requires calculating chi^2 gradient.
"""

import sys
import os
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl

from MulensModel import Event, Model, MulensData, MODULE_PATH


def chi2_fun(theta, event, parameters_to_fit):
    """
    for given event set attributes from parameters_to_fit (list of
    str) to values from theta list
    """
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return event.get_chi2()

def jacobian(theta, event, parameters_to_fit):
    """
    Calculate chi^2 gradient (also called Jacobian).
    """
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return event.chi2_gradient(parameters_to_fit)

# Read in the data file
data = MulensData(file_name=os.path.join(MODULE_PATH, "data", 
        "photometry_files", "phot_ob08092_O4.dat"))

# Initialize the fit
parameters_to_fit = ["t_0", "u_0", "t_E"]
t_0 = 5380.
u_0 = 0.1
t_E = 18.
model = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Link the data and the model
ev = Event(datasets=data, model=model)
print('Initial Trial\n{0}'.format(ev.model.parameters))

# Find the best-fit parameters
initial_guess = [t_0, u_0, t_E]
result = op.minimize(
    chi2_fun, x0=initial_guess, args=(ev, parameters_to_fit),
    method='Newton-CG',
    jac=jacobian, tol=1e-3)
(fit_t_0, fit_u_0, fit_t_E) = result.x

# Save the best-fit parameters
chi2 = chi2_fun(result.x, ev, parameters_to_fit)

# Output the fit parameters
msg = 'Best Fit: t_0 = {0:12.5f}, u_0 = {1:6.4f}, t_E = {2:8.3f}'
print(msg.format(fit_t_0, fit_u_0, fit_t_E))
print('Chi2 = {0:12.2f}'.format(chi2))
print('\nscipy.optimize.minimize result:')
print(result)

# Plot and compare the two models
init_model = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
final_model = Model({'t_0': fit_t_0, 'u_0': fit_u_0, 't_E': fit_t_E})
pl.figure()
init_model.plot_lc(data_ref=data, label='Initial Trial')
final_model.plot_lc(data_ref=data, label='Final Model')
pl.title('Difference b/w Input and Fitted Model')
pl.legend(loc='best')

# Plot the fitted model with the data
pl.figure()
ev.plot_data()
ev.plot_model(color='red')
pl.title('Data and Fitted Model')

pl.show()

