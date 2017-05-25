import scipy.optimize as op

import MulensModel

"""
Use Case 12: Fit a point lens to some data.
"""

#Initial Model
t_0 = 2457520.
u_0 = 0.6
t_E = 130.

model = MulensModel.Model()
model.set_parameters(t_0=t_0, u_0=u_0, t_E=t_E)

#Import data
data=MulensModel.MulensData(file_name='../../data/phot_ob160023.dat')

#Create Event
event = MulensModel.Event(datasets=data, model=model)
print('Initial Model')
print(event.model.parameters)

def lnlike(theta, event, parameters_to_fit):
    """
    likelihood function
    """
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model, val, theta[key])
    return event.get_chi2()

#Fit model to data using scipy
parameters_to_fit = ["t_0", "u_0", "t_E"]
result = op.minimize(lnlike, [t_0, u_0, t_E], args=(event, parameters_to_fit))
fit_t_0, fit_u_0, fit_t_E = result.x

#Save results
for key, val in enumerate(parameters_to_fit):
    setattr(event.model, val, result.x[key])

print('Fitted Model')
print(event.model.parameters)
