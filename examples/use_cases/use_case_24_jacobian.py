"""
Some scipy.optimize.minimize methods require a Jacobian to run. Here
is an example of how one would implement such a minimization method
using MulensModel.

Same as example_02_fitting.py except using the 'Newton-CG" method to
minimize the function (and now has a "Minimizer" class).

"""
import sys
import os
import numpy as np
import scipy.optimize as op
import matplotlib.pyplot as pl

import MulensModel
from MulensModel import Event, Fit, Model, MulensData, Utils


class Minimizer():
    """
    An object to link an Event to the functions necessary to minimize chi2.
    """

    def __init__(self, event, parameters_to_fit):
        self.event = event
        self.parameters_to_fit = parameters_to_fit

    def set_parameters(self, theta):
        """for given event set attributes from parameters_to_fit (list of str) 
        to values from theta list"""
        for (key, val) in enumerate(self.parameters_to_fit):
            setattr(self.event.model, val, theta[key])

    def chi2_fun(self, theta):
        """for a given set of parameters (theta), return the chi2"""
        self.set_parameters(theta)
        return self.event.get_chi2()
    
    def jacobian(self, theta):
        """for a given set of parameters (theta), return the jacobian"""
        self.set_parameters(theta) #might be redundant, but probably safer
        return self.event.jacobian()

#Read in the data file
SAMPLE_FILE_01 = os.path.join(MulensModel.MODULE_PATH, "data", 
                                                        "phot_ob08092_O4.dat")
data = MulensData(file_name=SAMPLE_FILE_01)

#Initialize the fit
parameters_to_fit = ["t_0", "u_0", "t_E"]
t_0 = 5380.
u_0 = 0.5
t_E = 18.
model = Model(t_0=t_0, u_0=u_0, t_E=t_E)

#Link the data and the model
ev = Event(datasets=data, model=model)
print('Initial Trial\n{0}'.format(ev.model.parameters))

#Create an object to hold the various minimization routines
minimizer = Minimizer(ev, parameters_to_fit)

#Find the best-fit parameters
initial_guess = [t_0, u_0, t_E]
result = op.minimize(
    minimizer.chi2_fun, x0=initial_guess, method='Newton-CG', 
    jac=minimizer.jacobian)

print(result.x)
(fit_t_0, fit_u_0, fit_t_E) = result.x

#Save the best-fit parameters
chi2 = minimizer.chi2_fun(result.x)

#Output the fit parameters
msg = 'Best Fit: t_0 = {0:12.5f}, u_0 = {1:6.4f}, t_E = {2:8.3f}'
print(msg.format(fit_t_0, fit_u_0, fit_t_E))
print('Chi2 = {0:12.2f}'.format(chi2))
print('scipy.optimize.minimize result:')
print(result)

#Plot and compare the two models
init_model = Model(t_0=t_0, u_0=u_0, t_E=t_E)
final_model = Model(t_0=fit_t_0, u_0=fit_u_0, t_E=fit_t_E)
pl.figure()
init_model.plot_lc(data_ref=data, label='Initial Trial')
final_model.plot_lc(data_ref=data, label='Final Model')
pl.title('Difference b/w Input and Fitted Model')
pl.legend(loc='best')

#Plot the fitted model with the data
pl.figure()
ev.plot_model()
ev.plot_data()
pl.title('Data and Fitted Model')

pl.show()

