"""
Use Case 15: Fit a point lens event with emcee.
"""
import numpy as np
import emcee
import os

import MulensModel as mm


def lnlike(theta, event, parameters_to_fit):
    """likelihood function """
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return -0.5 * (event.get_chi2() - chi2_0)


def lnprior(theta, parameters_to_fit):
    """priors"""
    if theta[parameters_to_fit.index("t_E")] < 0.:
        return -np.inf
    return 0.0


def lnprob(theta, event, parameters_to_fit):
    """combines likelihood and priors"""
    lp = lnprior(theta, parameters_to_fit)
    if not np.isfinite(lp):
        return -np.inf
    ln_like = lnlike(theta, event, parameters_to_fit)
    if np.isnan(ln_like):  # In the cases that source fluxes are negative we
        return -np.inf  # want to return these as if they were not in priors.
    return lp + ln_like


# Initialize the model
parameters_to_fit = ["t_0", "u_0", "t_E"]
parameters_values = [2455400., 0.5, 30.]
parameters_steps = [1., 0.01, 1.]

model = mm.Model(
    {'t_0': parameters_values[0], 'u_0': parameters_values[1],
     't_E': parameters_values[2]})
print("Initial", model.parameters)

# Read in the data
file_ = os.path.join(mm.DATA_PATH, "photometry_files",
                     "OB08092", "phot_ob08092_O4.dat")
data = mm.MulensData(file_name=file_, add_2450000=True)

# Set up the Event
event = mm.Event(datasets=data, model=model)

# Baseline chi2 = # of data points
chi2_0 = len(data.time) * 1.

# Initializations for emcee
ndim = len(parameters_values)
nwalkers = 100
nsteps = 500
burn = 50

start = [parameters_values + np.random.randn(ndim) * parameters_steps
         for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, lnprob, args=(event, parameters_to_fit))

# Run emcee - Fails because tries to set negative t_E.
# Verbose option to diagnose?
sampler.run_mcmc(start, nsteps)

samples = sampler.chain[:, burn:, :].reshape((-1, ndim))

results = map(
    lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
    zip(*np.percentile(samples, [16, 50, 84], axis=0)))

# Output fit
for r in results:
    print(*r)
