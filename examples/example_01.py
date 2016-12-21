#! /usr/bin/env python
"""
Example usage of MulensModel to fit a point lens light curve to the
data file phot_ob08092_04.dat.
"""
import sys
import numpy as np
import scipy.optimize as op

from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit
from MulensModel.event import Event
from MulensModel.model import Model
from MulensModel.utils import Utils


#Read in the data file
for path in sys.path:
    if path.find("MulensModel/source") > 0:
        MODULE_PATH = "/".join(path.split("/source")[:-1])
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"
data = MulensData(
    file_name=SAMPLE_FILE_01, date_fmt='HJDprime')

#Initialize the fit
parameters_to_fit = ["t_0", "u_0", "t_E"]
t_0 = 5380.
u_0 = 0.5
t_E = 18.
model = Model(t_0=t_0, u_0=u_0, t_E=t_E, coords="17:47:29.42 -34:43:35.6")

#Link the data and the model
ev = Event(datasets=data, model=model)

#Define the likelihood function
def lnlike(theta, event, parameters_to_fit):
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model, val, theta[key])
    return event.get_chi2()
        
#Find the best-fit parameters
result = op.minimize(lnlike, x0=[t_0, u_0, t_E], 
        args=(ev, parameters_to_fit), method='Nelder-Mead')
fit_t_0, fit_u_0, fit_t_E = result.x

#Save the best-fit parameters
for key, val in enumerate(parameters_to_fit):
    setattr(ev.model, val, result.x[key])
chi2 = ev.get_chi2()

#Output the fit parameters
print(
    'Best Fit: t_0 = {0:12.5f}, u_0 = {1:6.4f}, t_E = {2:8.3f}'.format(
        fit_t_0, fit_u_0, fit_tE))
print('Chi2 = {0:12.2f}'.format(chi2))
print('op.minimize result:')
print(result)

