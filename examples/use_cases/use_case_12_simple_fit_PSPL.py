"""
Use Case 12: Fit a point lens to some data.
"""
import os
import scipy.optimize as op

import MulensModel


#Initial Model
t_0 = 2457520.
u_0 = 0.6
t_E = 130.

model = MulensModel.Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

#Import data
file_name = os.path.join(
    MulensModel.MODULE_PATH, 'data', 'photometry_files', 'phot_ob160023.dat')
data = MulensModel.MulensData(file_name=file_name)

#Create Event
event = MulensModel.Event(datasets=data, model=model)
print('Initial Model')
print(event.model.parameters)

def chi2(theta, event, parameters_to_fit):
    """for given event set attributes from parameters_to_fit (list of str) 
    to values from theta list"""
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return event.get_chi2()

#Fit model to data using scipy
parameters_to_fit = ["t_0", "u_0", "t_E"]
initial_guess = [t_0, u_0, t_E]
result = op.minimize(chi2, initial_guess, args=(event, parameters_to_fit))
(fit_t_0, fit_u_0, fit_t_E) = result.x

#Save results and print.
chi2(result.x, event, parameters_to_fit)
print('Fitted Model')
print(event.model.parameters)
