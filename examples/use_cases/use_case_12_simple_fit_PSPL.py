"""
Use Case 12: Fit a point lens to some data.
"""
import os
import scipy.optimize as op

import MulensModel as mm


# Initial Model
t_0 = 2455380.
u_0 = 0.523
t_E = 20.

model = mm.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Import data
file_name = os.path.join(mm.DATA_PATH, 'photometry_files',
                         'OB08092', 'phot_ob08092_O4.dat')
data = mm.MulensData(file_name=file_name, add_2450000=True)

# Create Event
event = mm.Event(datasets=data, model=model)
print('Initial Model')
print(event.model.parameters)
print(event.get_chi2())


def chi2(theta, event, parameters_to_fit):
    """
    for given event set attributes from parameters_to_fit (list of str)
    to values from theta list
    """
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])

    return event.get_chi2()


# Fit model to data using scipy
parameters_to_fit = ["t_0", "u_0", "t_E"]
initial_guess = [t_0, u_0, t_E]
result = op.minimize(
    chi2, initial_guess, args=(event, parameters_to_fit),
    method='Nelder-Mead')
(fit_t_0, fit_u_0, fit_t_E) = result.x

# Save results and print.
chi2(result.x, event, parameters_to_fit)
print('Fitted Model')
print(event.model.parameters)
print(event.get_chi2())
