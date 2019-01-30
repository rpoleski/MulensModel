import os
import sys
import numpy as np
import emcee
from matplotlib import pyplot as plt
import configparser

import MulensModel as MM

import hack_read as read


def ln_like(theta, event, parameters_to_fit, print_models):
    """
    Likelihood function. The values of *parameters_to_fit* are in *theta*.
    MulensModel Event class instance *event* gives event for which
    calculations will be done. Boolean *print_models* controls if
    all models are printed.
    """
    for (theta_, param) in zip(theta, parameters_to_fit):
        setattr(event.model.parameters, param, theta_)
    chi2 = event.get_chi2()
    if print_models:
        print(chi2, *[t for t in theta], flush=True)
    return -0.5 * chi2

def ln_prior(theta, parameters_to_fit):
    """
    Prior. Check if *theta* values for *parameters_to_fit* are within ranges
    defined by *ln_prior.min* and *ln_prior.max*.
    """
    inside = 0.
    outside = -np.inf

    for (parameter, value) in ln_prior.min.items():
        index = parameters_to_fit.index(parameter)
        if theta[index] < value:
            return outside

    for (parameter, value) in ln_prior.max.items():
        index = parameters_to_fit.index(parameter)
        if theta[index] > value:
            return outside

    return inside

def ln_prob(
        theta, event, parameters_to_fit, print_models=False):
    """
    Log probability of the model - combines ln_prior() and ln_like().
    """
    ln_prior_ = ln_prior(theta, parameters_to_fit)
    if not np.isfinite(ln_prior_):
        return -np.inf

    ln_like_ = ln_like(theta, event, parameters_to_fit, print_models)
    if np.isnan(ln_like_):
        return -np.inf

    return ln_prior_ + ln_like_

def generate_random_parameters(parameters, starting, n):
    """
    Generate *n* vectors of values of *parameters* according to distributions
    specified in *starting*.
    """
    values = []
    for param in parameters:
        settings = starting[param]
        if settings[0] == 'gauss':
            v = settings[2] * np.random.randn(n)
            v += settings[1]
        elif settings[0] == 'uniform':
            v = np.random.uniform(
                low=settings[1], high=settings[2], size=n)
        elif settings[0] == 'log-uniform':
            beg = np.log(settings[1])
            end = np.log(settings[2])
            v = np.exp(np.random.uniform(beg, end, n))
        values.append(v)
    return np.array(values).T.tolist()


# Read config file.
if len(sys.argv) != 2:
    raise ValueError('Exactly one argument needed - cfg file')
config_file = sys.argv[1]

config = configparser.ConfigParser()
config.optionxform = str # So that "t_E" is not changed to "t_e".
config.read(config_file)
files = read.read_files_from_config(config)
model_settings = read.read_model_settings(config)
(parameters, starting) = read.read_parameters_start(config)
fixed_parameters = read.read_fix_parameters(config)
(min_values, max_values) = read.read_min_max(config)
ln_prior.min = min_values
ln_prior.max = max_values
emcee_settings = read.read_emcee_settings(config)
other_settings = read.read_other(config)

# Read photometric data.
datasets = [MM.MulensData(file_name=f[0], phot_fmt=f[1]) for f in files]

# Generate starting values of parameters.
start = generate_random_parameters(parameters, starting, emcee_settings['n_walkers'])

# Setup Event instance that combines model and data.
par = dict(zip(parameters, start[0]))
par = {**par, **fixed_parameters}
my_model = MM.Model(par, coords=model_settings['coords'])
if 'methods' in model_settings:
    my_model.set_magnification_methods(model_settings['methods'])
if 'default_method' in model_settings:
    my_model.set_default_magnification_method(model_settings['default_method'])
my_event = MM.Event(datasets=datasets, model=my_model)

# Prepare sampler.
n_dim = len(parameters)
print_models = other_settings.get('print_models', False)
args = (my_event, parameters, print_models)
sampler = emcee.EnsembleSampler(emcee_settings['n_walkers'], n_dim, ln_prob, args=args)

# Run sampler.
sampler.run_mcmc(start, emcee_settings['n_steps'])

# Parse results.
burn = emcee_settings['n_burn']
samples = sampler.chain[:, burn:, :].reshape((-1, n_dim))
r_16 = np.percentile(samples, 16, axis=0)
r_50 = np.percentile(samples, 50, axis=0)
r_84 = np.percentile(samples, 84, axis=0)
print("Fitted parameters:")
for i in range(n_dim):
    if parameters[i] == 'q':
        fmt = "{:} {:.7f} +{:.7f} -{:.7f}"
    else:
        fmt = "{:} {:.5f} +{:.5f} -{:.5f}"
    print(fmt.format(parameters[i], r_50[i], r_84[i]-r_50[i], r_50[i]-r_16[i]))
print("Smallest chi2 model:")
best = [my_event.best_chi2_parameters[p] for p in parameters]
print(*[b if isinstance(b, float) else b.value for b in best])
print(my_event.best_chi2)

# Plot results.
ln_like(best, my_event, parameters, False) # This allows plotting of the best model.
print(my_event.model)
my_event.plot_data(subtract_2450000=True)
my_event.plot_model(
    subtract_2450000=True,
    t_start=other_settings['plot_time'][0]+2450000.,
    t_stop=other_settings['plot_time'][1]+2450000.)
plt.xlim(*other_settings['plot_time'])
plt.show()
