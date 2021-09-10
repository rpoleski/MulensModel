"""
Fit point source point lens model to OB08092 using Newton-CG method from scipy.
This method requires calculating chi^2 gradient.
"""

import sys
import os
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as plt

import MulensModel as mm


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

    Note: this implementation is robust but possibly inefficient. If
    chi2_fun() is ALWAYS called before jacobian with the same parameters,
    there is no need to set the parameters in event.model; also,
    event.calculate_chi2_gradient() can be used instead (which avoids fitting
    for the fluxes twice).
    """
    for (key, val) in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])

    return event.get_chi2_gradient(parameters_to_fit)


# Read in the data file
file_ = os.path.join(mm.DATA_PATH, "photometry_files", "OB08092",
                     "phot_ob08092_O4.dat")
data = mm.MulensData(file_name=file_)

# Initialize the fit
parameters_to_fit = ["t_0", "u_0", "t_E"]
t_0 = 5380.
u_0 = 0.1
t_E = 18.
model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Link the data and the model
ev = mm.Event(datasets=data, model=model)
(source_flux_init, blend_flux_init) = ev.get_flux_for_dataset(data)
print('Initial Trial\n{0}'.format(ev.model.parameters))
print('Chi2 = {0}\n'.format(ev.get_chi2()))

# Find the best-fit parameters
initial_guess = [t_0, u_0, t_E]
result = op.minimize(
    chi2_fun, x0=initial_guess, args=(ev, parameters_to_fit),
    method='Newton-CG',
    jac=jacobian, tol=1e-3)
(fit_t_0, fit_u_0, fit_t_E) = result.x

# Save the best-fit parameters
chi2 = chi2_fun(result.x, ev, parameters_to_fit)
(source_flux_final, blend_flux_final) = ev.get_flux_for_dataset(data)

# Output the fit parameters
msg = 'Best Fit: t_0 = {0:12.5f}, u_0 = {1:6.4f}, t_E = {2:8.3f}'
print(msg.format(fit_t_0, fit_u_0, fit_t_E))
print('Chi2 = {0:12.2f}'.format(chi2))
print('\nscipy.optimize.minimize result:')
print(result)

# Plot and compare the two models
init_model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})
final_model = mm.Model({'t_0': fit_t_0, 'u_0': fit_u_0, 't_E': fit_t_E})
plt.figure()
init_model.plot_lc(
    source_flux=source_flux_init, blend_flux=blend_flux_init,
    label='Initial Trial')
final_model.plot_lc(
    source_flux=source_flux_final, blend_flux=blend_flux_final,
    label='Final Model')
plt.title('Difference b/w Input and Fitted Model')
plt.legend(loc='best')

# Plot the fitted model with the data
plt.figure()
ev.plot_data()
ev.plot_model(color='red')
plt.title('Data and Fitted Model')

plt.show()
