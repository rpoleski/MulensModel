"""
Fits PSPL model with parallax ground-based and satellite data
using EMCEE sampler. We're using photometry of OGLE-2014-BLG-0939 from:

Yee et al. 2015 ApJ 802, 76
https://ui.adsabs.harvard.edu/abs/2015ApJ...802...76Y/abstract

and explore only 1 of 4 degenerate models.

It is similar to example_06_fit_parallax_EMCEE.py
"""
from os.path import join as join
import sys
import numpy as np
try:
    import emcee
except ImportError as err:
    print(err)
    print("\nEMCEE could not be imported.")
    print("Get it from: http://dfm.io/emcee/current/user/install/")
    print("and re-run the script")
    sys.exit(1)
import matplotlib.pyplot as plt

import MulensModel as mm


# Define likelihood functions
def ln_like(theta, event, parameters_to_fit):
    """likelihood function"""
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return -0.5 * event.get_chi2()


def ln_prior(theta, parameters_to_fit):
    """priors - we only reject obviously wrong models"""
    if theta[parameters_to_fit.index("t_E")] < 0.:
        return -np.inf
    return 0.0


def ln_prob(theta, event, parameters_to_fit):
    """combines likelihood and priors"""
    ln_prior_ = ln_prior(theta, parameters_to_fit)
    if not np.isfinite(ln_prior_):
        return -np.inf
    ln_like_ = ln_like(theta, event, parameters_to_fit)

    # In the cases that source fluxes are negative we want to return
    # these as if they were not in priors.
    if np.isnan(ln_like_):
        return -np.inf

    return ln_prior_ + ln_like_


# Read the data (note that we do not rescale errorbars here):
dir_ = join(mm.DATA_PATH, "photometry_files", "OB140939")
file_ground = join(dir_, "ob140939_OGLE.dat")
file_spitzer = join(dir_, "ob140939_Spitzer.dat")
data_ground = mm.MulensData(
    file_name=file_ground, plot_properties={'label': 'OGLE'})

# Here is the main difference - we provide the ephemeris for Spitzer:
file_spitzer_eph = join(
    mm.DATA_PATH, 'ephemeris_files', 'Spitzer_ephemeris_01.dat')
data_spitzer = mm.MulensData(
    file_name=file_spitzer, ephemerides_file=file_spitzer_eph,
    plot_properties={'label': 'Spitzer'})

# For parallax calculations we need event coordinates:
coords = "17:47:12.25 -21:22:58.7"

# Starting parameters:
params = {
    't_0': 2456830., 'u_0': 0.8, 't_E': 25.,
    'pi_E_N': 0., 'pi_E_E': 0.,
    't_0_par': 2456836.06}
my_model = mm.Model(params, coords=coords)
my_event = mm.Event(datasets=[data_ground, data_spitzer], model=my_model)

# Which parameters we want to fit?
parameters_to_fit = ["t_0", "u_0", "t_E", "pi_E_N", "pi_E_E"]
# And remember to provide dispersions to draw starting set of points
sigmas = [0.1, 0.01, 0.1, 0.05, 0.05]

# Initializations for EMCEE
n_dim = len(parameters_to_fit)
n_walkers = 20
n_steps = 1500
n_burn = 500
# Including the set of n_walkers starting points:
start_1 = [params[p] for p in parameters_to_fit]
start = [start_1 + np.random.randn(n_dim) * sigmas
         for i in range(n_walkers)]

# Run emcee (this should take about a minute):
sampler = emcee.EnsembleSampler(
    n_walkers, n_dim, ln_prob, args=(my_event, parameters_to_fit))
sampler.run_mcmc(start, n_steps)

# Remove burn-in samples and reshape:
samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))

# Results:
results = np.percentile(samples, [16, 50, 84], axis=0)
print("Fitted parameters:")
fmt = "{:} : {:.5f} {:.5f} {:.5f}"
for (i, p) in enumerate(parameters_to_fit):
    r = results[1, i]
    print(fmt.format(p, r, results[2, i]-r, r-results[0, i]))

# We extract best model parameters and chi2 from the chain:
prob = sampler.lnprobability[:, n_burn:].reshape((-1))
best_index = np.argmax(prob)
best_chi2 = prob[best_index] / -0.5
best = samples[best_index, :]
print("\nSmallest chi2 model:")
print(*[repr(b) if isinstance(b, float) else b.value for b in best])
print(best_chi2)
for (i, parameter) in enumerate(parameters_to_fit):
    setattr(my_event.model.parameters, parameter, best[i])

my_event.fit_fluxes()

# In order to make plots, we need a Model instance
# that has satellite ephemeris:
params = my_event.model.parameters.parameters
space_model = mm.Model({**params}, coords=coords,
                       ephemerides_file=file_spitzer_eph)

# Prepare plots:
my_event.plot_model(subtract_2450000=True)
fluxes = my_event.get_flux_for_dataset(data_ground)
# We need this to ensure that fluxes are scaled properly.
space_model.plot_lc(subtract_2450000=True,
                    source_flux=fluxes[0], blend_flux=fluxes[1])
my_event.plot_data(subtract_2450000=True)
plt.legend()
plt.xlim(6800., 6880.)

plt.figure()
my_event.model.plot_trajectory()
space_model.plot_trajectory()
space_model.plot_caustics(color='black')
plt.axis('equal')
plt.xlim(-1.1, 1.1)
plt.show()
