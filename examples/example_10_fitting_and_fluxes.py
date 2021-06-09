"""
Fits PSPL model for MB08310 with finite source effect using EMCEE sampler and
reports posterior flux parameters. The flux parameters are found by regression
for each model and cached using EMCEE blobs mechanism.

This example uses config file:

python example_10_fitting_and_fluxes.py example_10.cfg

Data from `Janczak et al. 2010, ApJ 711, 731
<https://ui.adsabs.harvard.edu/abs/2010ApJ...711..731J/abstract>`_.

"""

import sys
from os.path import isfile
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
import configparser

import MulensModel as mm


# Define likelihood functions
def ln_like(theta, event, parameters_to_fit):
    """ likelihood function """
    for key, val in enumerate(parameters_to_fit):
        setattr(event.model.parameters, val, theta[key])
    return -0.5 * event.get_chi2()


def ln_prior(theta, parameters_to_fit):
    """priors - we only reject obviously wrong models"""
    if theta[parameters_to_fit.index("t_E")] < 0.:
        return -np.inf
    if theta[parameters_to_fit.index("t_star")] < 0.:
        return -np.inf
    return 0.0


def get_fluxes(event):
    """
    For given Event instance extract all the fluxes and return them in
    a list. Odd elements are source fluxes and even ones are blending fluxes.
    These fluxes are in units used by MulensModel, where 1 corresponds to
    22 mag.
    """
    fluxes = []
    for dataset in event.datasets:
        (data_source_flux, data_blend_flux) = event.get_flux_for_dataset(
            dataset)
        fluxes.append(data_source_flux[0])
        fluxes.append(data_blend_flux)

    return fluxes


def ln_prob(theta, event, parameters_to_fit):
    """
    combines likelihood and priors; returns ln(prob) and a list of fluxes
    """
    fluxes = [None] * 2 * len(event.datasets)

    ln_prior_ = ln_prior(theta, parameters_to_fit)
    if not np.isfinite(ln_prior_):
        return (-np.inf, fluxes)
    ln_like_ = ln_like(theta, event, parameters_to_fit)

    # In the cases that source fluxes are negative we want to return
    # these as if they were not in priors.
    if np.isnan(ln_like_):
        return (-np.inf, fluxes)

    ln_prob_ = ln_prior_ + ln_like_
    fluxes = get_fluxes(event)
    return (ln_prob_, fluxes)


def uncertainties(x):
    """
    change the parameters so that we can easily print median and uncertainties
    (based on percentiles)
    """
    return (x[1], x[2]-x[1], x[1]-x[0])


# Open config file
if len(sys.argv) != 2:
    raise ValueError('Exactly one argument needed - cfg file')
config_file = sys.argv[1]
if not isfile(config_file):
    raise FileNotFoundError('File: {:}'.format(config_file))
config = configparser.ConfigParser()
config.optionxform = str
config.read(config_file)

# Read the data
section = "photometry files"
if section not in config:
    raise KeyError('Sorry, no photometry files specified in config.')
file_names = [config.get(section, var) for var in config[section]]
kwargs = {'comments': ["\\", "|"]}
data = [mm.MulensData(file_name=name, **kwargs) for name in file_names]

# Read parameters
section = "parameters to fit"
info = [[var, config.get(section, var).split()] for var in config[section]]
for info_ in info:
    if len(info_[1]) != 2:
        msg = 'Wrong input in cfg file:\n{:}'
        raise ValueError(mag.format(config.get(section, info_[0])))
parameters_to_fit = [x[0] for x in info]
starting_mean = [float(x[1][0]) for x in info]
starting_sigma = [float(x[1][1]) for x in info]
parameters = {key: va for (key, va) in zip(parameters_to_fit, starting_mean)}
n_parameters = len(parameters)

model = mm.Model(parameters)
event = mm.Event(datasets=data, model=model)

# Read methods
section = "methods"
if section in config:
    text = config.get(section, section).split()
    methods = [t if i % 2 == 1 else float(t) for (i, t) in enumerate(text)]
    model.set_magnification_methods(methods)

# Read EMCEE options
section = "EMCEE parameters"
n_walkers = config.getint(section, 'walkers')
n_steps = config.getint(section, 'steps')
if 'burn' in config[section]:
    n_burn = config.getint(section, 'burn')
else:
    n_burn = 0

# Read plotting kwargs
section = "plot kwargs"
plot_kwargs = dict()
if section in config:
    for key in config[section]:
        plot_kwargs[key] = config.getfloat(section, key)

# Starting points:
starting = [starting_mean + starting_sigma * np.random.randn(n_parameters)
            for i in range(n_walkers)]

# Run the sampler
sampler = emcee.EnsembleSampler(
    n_walkers, n_parameters, ln_prob, args=(event, parameters_to_fit))
sampler.run_mcmc(starting, n_steps)

samples = sampler.chain[:, n_burn:, :].reshape((-1, n_parameters))
# The sampler.blobs is where EMCEE passes all the additional parameters
# returned by ln_prob().
blob_sampler = np.transpose(np.array(sampler.blobs), axes=(1, 0, 2))
n_fluxes = blob_sampler.shape[-1]
blob_samples = blob_sampler[:, n_burn:, :].reshape((-1, n_fluxes))

# Results:
results = np.percentile(samples, [16, 50, 84], axis=0)
print("Fitted parameters:")
for i in range(n_parameters):
    print("{:.8s} : {:.5f} {:.5f} {:.5f}".format(parameters_to_fit[i],
          *uncertainties(results[:, i])))
blob_results = np.percentile(blob_samples, [16, 50, 84], axis=0)
flux_name = ['s', 'b']
print("Fluxes (source and blending):")
for i in range(n_fluxes):
    print('flux_{:}_{:} : {:.4f} {:.4f} {:.4f}'.format(flux_name[i % 2],
          i//2+1, *uncertainties(blob_results[:, i])))

# We extract best model parameters and chi2:
print("\nSmallest chi2 model:")
prob = sampler.lnprobability[:, n_burn:].reshape((-1))
best_index = np.argmax(prob)
best_chi2 = prob[best_index] / -0.5
best = samples[best_index, :]
for key, val in enumerate(parameters_to_fit):
    setattr(event.model.parameters, val, best[key])

print("\nSmallest chi2 model:")
print(*[repr(b) if isinstance(b, float) else b.value for b in best])
print("chi2 = ", event.get_chi2())
# Note: the call to event.get_chi2() is necessary so that event is updated with
# the fluxes for the best-fit model for use in plotting (below).

# Plot model and data.
print("\nNow let's plot the best model")
event.plot_model(subtract_2450000=True, **plot_kwargs)
ylim = plt.ylim()
# The command below raises warning - it's not your fault. It's caused
# by the input data.
event.plot_data(subtract_2450000=True)
xlim = list(plt.xlim())
if 't_start' in plot_kwargs:
    xlim[0] = plot_kwargs['t_start'] - 2450000.
if 't_stop' in plot_kwargs:
    xlim[1] = plot_kwargs['t_stop'] - 2450000.
plt.ylim(*ylim)
plt.xlim(*xlim)
plt.legend(loc='best')
plt.show()
