"""
Fits binary source model using EMCEE sampler.

The code simulates binary source light curve and fits the model twice:
with source flux ratio found via linear regression and
with source flux ratio as a chain parameter.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
try:
    import emcee
except ImportError as err:
    print(err)
    print("\nEMCEE could not be imported.")
    print("Get it from: http://dfm.io/emcee/current/user/install/")
    print("and re-run the script")
    sys.exit(1)

import MulensModel as mm


# Fix the seed for the random number generator so the behavior is reproducible.
np.random.seed(12343)


# Define likelihood functions
def ln_like(theta, event, parameters_to_fit):
    """ likelihood function """
    for (param, theta_) in zip(parameters_to_fit, theta):
        # Here we handle fixing source flux ratio:
        if param == 'flux_ratio':
            # implemented for a single dataset
            event.fix_source_flux_ratio = {my_dataset: theta_}
        else:
            setattr(event.model.parameters, param, theta_)

    return -0.5 * event.get_chi2()


def ln_prior(theta, parameters_to_fit):
    """priors - we only reject obviously wrong models"""
    for param in ['t_E', 'u_0_1', 'u_0_2']:
        if param in parameters_to_fit:
            if theta[parameters_to_fit.index(param)] < 0.:
                return -np.inf
    return 0.0


def ln_prob(theta, event, parameters_to_fit):
    """ combines likelihood and priors"""
    ln_prior_ = ln_prior(theta, parameters_to_fit)
    if not np.isfinite(ln_prior_):
        return -np.inf
    ln_like_ = ln_like(theta, event, parameters_to_fit)

    # In the cases that source fluxes are negative we want to return
    # these as if they were not in priors.
    if np.isnan(ln_like_):
        return -np.inf

    return ln_prior_ + ln_like_


def fit_EMCEE(parameters_to_fit, starting_params, sigmas, ln_prob, event,
              n_walkers=20, n_steps=3000, n_burn=1500):
    """
    Fit model using EMCEE and print results.
    Arguments:
        parameters_to_fit - list of parameters
        starting_params - dict that specifies values of these parameters
        sigmas - list of sigma values used to find starting values
        ln_prob - function returning logarithm of probability
        event - MulensModel.Event instance
        n_walkers - number of walkers in EMCEE
        n_steps - number of steps per walker
        n_burn - number of steps considered as burn-in ( < n_steps)
    """
    n_dim = len(parameters_to_fit)
    mean = [starting_params[p] for p in parameters_to_fit]
    start = [mean + np.random.randn(n_dim) * sigmas for i in range(n_walkers)]

    # Run emcee (this can take some time):
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, ln_prob, args=(event, parameters_to_fit))
    sampler.run_mcmc(start, n_steps)

    # Remove burn-in samples and reshape:
    samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))

    # Results:
    results = np.percentile(samples, [16, 50, 84], axis=0)
    print("Fitted parameters:")
    for i in range(n_dim):
        r = results[1, i]
        msg = parameters_to_fit[i] + ": {:.5f} +{:.5f} -{:.5f}"
        print(msg.format(r, results[2, i]-r, r-results[0, i]))

    # We extract best model parameters and chi2 from event:
    prob = sampler.lnprobability[:, n_burn:].reshape((-1))
    best_index = np.argmax(prob)
    best = samples[best_index, :]
    for key, val in enumerate(parameters_to_fit):
        if val == 'flux_ratio':
            event.fix_source_flux_ratio = {my_dataset: best[key]}
        else:
            setattr(event.model.parameters, val, best[key])

    print("\nSmallest chi2 model:")
    print(*[repr(b) if isinstance(b, float) else b.value for b in best])
    print("chi2 = ", event.get_chi2())


# First, prepare the data. There is nothing very exciting in this part,
# so you may skip it.
t_0_1 = 6100.
u_0_1 = 0.2
t_0_2 = 6140.
u_0_2 = 0.01
t_E = 25.
assumed_flux_1 = 100.
assumed_flux_2 = 5.
assumed_flux_blend = 10.
n_a = 1000
n_b = 600
time_a = np.linspace(6000., 6300., n_a)
time_b = np.linspace(6139., 6141., n_b)
time = np.sort(np.concatenate((time_a, time_b)))
model_1 = mm.Model({'t_0': t_0_1, 'u_0': u_0_1, 't_E': t_E})
A_1 = model_1.get_magnification(time)
model_2 = mm.Model({'t_0': t_0_2, 'u_0': u_0_2, 't_E': t_E})
A_2 = model_2.get_magnification(time)
flux = A_1 * assumed_flux_1 + A_2 * assumed_flux_2 + assumed_flux_blend
flux_err = 6. + 0. * time
flux += flux_err * np.random.normal(size=n_a+n_b)
my_dataset = mm.MulensData([time, flux, flux_err], phot_fmt='flux')
# If you want to plot, then just uncomment:
# plt.plot(time, flux, 'ro')
# plt.show()

# Starting parameters:
params = {'t_0_1': 6101., 'u_0_1': 0.19, 't_0_2': 6140.123, 'u_0_2': 0.04,
          't_E': 20.}
my_model = mm.Model(params)
my_event = mm.Event(datasets=my_dataset, model=my_model)

# First fit - source flux ratio not set, hence found by regression:
parameters_to_fit = ["t_0_1", "u_0_1", "t_0_2", "u_0_2", "t_E"]
sigmas = [0.1, 0.05, 1., 0.01, 1.]
print("\nFirst fit. This can take some time...")
fit_EMCEE(parameters_to_fit, params, sigmas, ln_prob, my_event)

# Starting parameters for second fit:
params = {'t_0_1': 6101., 'u_0_1': 0.19, 't_0_2': 6140.123, 'u_0_2': 0.04,
          't_E': 25.987}
my_model = mm.Model(params)
my_event = mm.Event(datasets=my_dataset, model=my_model)
params['flux_ratio'] = 0.02

# Second fit - source flux ratio as one of the chain parameters:
parameters_to_fit = ["t_0_1", "u_0_1", "t_0_2", "u_0_2", "t_E", "flux_ratio"]
sigmas = [0.1, 0.05, 1., 0.01, 1., 0.001]
print("\nSecond fit. This can take some time...")
fit_EMCEE(parameters_to_fit, params, sigmas, ln_prob, my_event)
