#! /usr/bin/env python

import numpy as np
import emcee

import MulensModel


def lnlike(theta, event, parameters_to_fit):
    for key, val in iterate(parameters_to_fit):
        setattr(event.model, val, theta[key])
    return -0.5 * (event.get_chi2() - event.chi2_0)

def lnprior(theta):
    return 0.0

def lnprob(theta, event, parameters_to_fit):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, event, parameters_to_fit)

parameters_to_fit = ["t_0", "u_0", "t_E"]
parameters_values = [7000.1, 0.1, 100.]

model = MulensModel.Model()
for key, val in iterate(parameters_to_fit):
    setattr(model, val, parameters_values[key])

data=MulensModel.MulensData(file_name="data_file.dat")

event = MulensModel.Event(datasets=data, model=model)
event.chi2_0 = len(data) * 1.

ndim = len(parameters_values)
nwalkers = 100
nsteps = 500
burn = 50

start = [parameters_values * (1 + 1e-4 * np.random.randn(ndim)) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(event, parameters_to_fit))

sampler.run_mcmc(start, nsteps)

samples = sampler.chain[:, burn:, :].reshape((-1, ndim))

results = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))

for r in results:
    print(*r)

