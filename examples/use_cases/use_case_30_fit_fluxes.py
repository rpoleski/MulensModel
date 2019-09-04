"""
Assume you want (some) of the fluxes to be freely fitted parameters, e.g.
so you can measure the uncertainty in the color.
"""
import os
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

import MulensModel

raise NotImplementedError('This use case has not been implemented.')

def ln_like(theta, event, parameters_to_fit):
    """ likelihood function """
    KCT01 = 0
    Spitzer = 9
    for (key, val) in enumerate(parameters_to_fit):
        # Some fluxes are MCMC parameters
        if val[0:] == 'f':
            if val == 'fsKCT01':
                event.datasets[KCT01].source_flux = theta[key]
            elif val == 'fsSpitzer':
                event.datasets[Spitzer].source_flux = theta[key]

        else:
            setattr(event.model.parameters, val, theta[key])

    return -0.5 * event.get_chi2()


def ln_prior(theta, parameters_to_fit):
    """priors - we only reject obviously wrong models"""
    if theta[parameters_to_fit.index("t_E")] < 0.:
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


# Add data from OB161195
datasets = []
file_names = ['KCT01I.dat', 'KCT41I.dat', 'KCT42I.dat', 'KSA01I.dat',
              'KSA41I.dat', 'KSA42I.dat', 'KSS01I.dat', 'KSS41I.dat',
              'KSS42I.dat', 'spitzer_b12.dat']
dir_ = os.path.join(MulensModel.DATA_PATH, "photometry_files", "OB161195")
for file_name in file_names:
    file_ = os.path.join(dir_, file_name)
    datasets.append(MulensModel.MulensData(file_name=file_, add_2450000=True))

# Close-- model
model = MulensModel.Model(
    {'t_0': 2457568.7692, 'u_0': -0.05321, 't_E': 9.96, 'rho': 0.00290,
     'pi_E_N': -0.2154, 'pi_E_E': -0.380,
     'alpha': np.rad2deg(-0.9684), 's': 0.9842, 'q': 0.0000543})

methods = [7560., 'VBBL', 7580.]

model.set_magnification_methods(methods)
model.set_default_magnification_method('point_source_point_lens')

event = MulensModel.Event(datasets=datasets, model=model,
                          coords="17:55:23.50 -30:12:26.1")

# Which parameters we want to fit?
parameters_to_fit = ["t_0", "u_0", "t_E", "pi_E_N", "pi_E_E", "fsKCT01", "fsSpitzer"]
# And remember to provide dispersions to draw starting set of points
sigmas = [0.01, 0.001, 0.1, 0.01, 0.01, 0.01, 0.01]

# Initializations for EMCEE
n_dim = len(parameters_to_fit)
n_walkers = 40
n_steps = 500
n_burn = 150
# Including the set of n_walkers starting points:
start_1 = [params[p] for p in parameters_to_fit]
start = [start_1 + np.random.randn(n_dim) * sigmas
         for i in range(n_walkers)]

# Run emcee (this can take some time):
sampler = emcee.EnsembleSampler(
    n_walkers, n_dim, ln_prob, args=(my_event, parameters_to_fit))
sampler.run_mcmc(start, n_steps)

# Remove burn-in samples and reshape:
samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))

# Results:
results = np.percentile(samples, [16, 50, 84], axis=0)
print("Fitted parameters:")
for i in range(n_dim):
    r = results[1, i]
    print("{:.5f} {:.5f} {:.5f}".format(r, results[2, i]-r, r-results[0, i]))

# We extract best model parameters and chi2 from my_event:
print("\nSmallest chi2 model:")
best = [my_event.best_chi2_parameters[p] for p in parameters_to_fit]
print(*[repr(b) if isinstance(b, float) else b.value for b in best])
print(my_event.best_chi2)

# Compare to the Color constraint for OB161195 (I_KMT - L_Spitzer)
(source_color, sigma_color) = (0.78, 0.03)
plt.hist(-2.5*np.log10(samples[:,0] / samples[:,9]), bins=25)
plt.axvline(source_color, linestyle='-', color='black')
plt.axvline(source_color - sigma_color, linestyle=':', color='black')
plt.axvline(source_color + sigma_color, linestyle=':', color='black')
pl.xlabel('(I-L)')
plt.show()

