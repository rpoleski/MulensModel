#! /usr/bin/env python

import numpy as np
import emcee

import MulensModel

"""
Use Case 15: Fit a point lens event with emcee. 
"""


def lnlike(theta, event, parameters_to_fit):
    """ likelihood function """
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model, val, theta[key])
    return -0.5 * (event.get_chi2() - chi2_0)

def lnprior(theta, parameters_to_fit):
    """priors"""
    if theta[parameters_to_fit.index("t_E")] < 0.:
        return -np.inf
    return 0.0

def lnprob(theta, event, parameters_to_fit):
    """ conbines likelihood and priors"""
    lp = lnprior(theta, parameters_to_fit)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, event, parameters_to_fit)

#Initialize the model
parameters_to_fit = ["t_0", "u_0", "t_E"]
parameters_values = [2457500., 0.5, 100.]

model = MulensModel.Model()
for key, val in enumerate(parameters_to_fit):
    setattr(model, val, parameters_values[key])
print("Initial", model.parameters)

#Read in the data
data = MulensModel.MulensData(file_name="../../data/phot_ob160023.dat")

#Set up the Event
event = MulensModel.Event(datasets=data, model=model)

#Baseline chi2 = # of data points
chi2_0 = len(data.time) * 1.

#Initializations for emcee
ndim = len(parameters_values)
nwalkers = 100
nsteps = 500
burn = 50

start = [parameters_values * (1. + 1e-4 * np.random.randn(ndim)) 
         for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, lnprob, args=(event, parameters_to_fit))

#Run emcee - Fails because tries to set negative tE. Verbose option to diagnose?
sampler.run_mcmc(start, nsteps)

samples = sampler.chain[:, burn:, :].reshape((-1, ndim))

results = map(
    lambda v: (v[1], v[2]-v[1], v[1]-v[0]), 
    zip(*np.percentile(samples, [16, 50, 84], axis=0)))

#Output fit
for r in results:
#    print(*r)
    print(r)

