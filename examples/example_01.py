#! /usr/bin/env python

import sys
import numpy as np
import scipy.optimize as op

from MulensModel.mulensdata import MulensData
from MulensModel.fit import Fit
from MulensModel.event import Event
from MulensModel.model import Model
from MulensModel.utils import Utils


for path in sys.path:
    if path.find("MulensModel/source") > 0:
        MODULE_PATH = "/".join(path.split("/source")[:-1])
SAMPLE_FILE_01 = MODULE_PATH + "/data/phot_ob08092_O4.dat"

parameters_to_fit = ["t_0", "u_0", "t_E"]
t_0 = 5380.
u_0 = 0.5
t_E = 18.
    
data = MulensData(file_name=SAMPLE_FILE_01)
    
mod = Model(t_0=t_0, u_0=u_0, t_E=t_E)
mod.set_datasets([data])
ev = Event()
#ev = Event(datasets=data, model=model)
ev.model = mod
ev.datasets = [data]


def lnlike(theta, event, parameters_to_fit):
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model, val, theta[key])
    return event.get_chi2()
        
result = op.minimize(lnlike, x0=[t_0, u_0, t_E], 
        args=(ev, parameters_to_fit), method='Nelder-Mead')
fit_t_0, fit_u_0, fit_t_E = result.x

for key, val in enumerate(parameters_to_fit):
    setattr(ev.model, val, result.x[key])
        
chi2 = ev.get_chi2()

print(chi2)
print(result["x"])
print(result)

