"""
This example fits model to OGLE-2003-BLG-235/MOA-2003-BLG-53,
the first microlensing planet. Here we fix *s* and *q* parameters for
the sake of simplicity. Wide range of other binary lens parameters is explored.

Note that it would be beneficial to turn *x_caustic_in* and *x_caustic_out*
to periodic variables.

Specific settings are in file example_13.cfg.

Running this example takes approx 5 minutes on most modern machines.
"""
import numpy as np
import emcee
import configparser

import MulensModel as mm

import example_15_read as read


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
        if 'x_caustic_in' not in parameters_to_fit:
            print(chi2, *[t for t in theta], flush=True)
        else:
            theta_ = theta.tolist()
            keys = ['t_0', 'u_0', 't_E', 'alpha']
            theta_ += [getattr(event.model.parameters, key) for key in keys]
            print(chi2, *theta_, flush=True)
    return -0.5 * chi2


def ln_prior(theta, parameters_to_fit, event):
    """
    Prior. Check if *theta* values for *parameters_to_fit* are within ranges
    defined by *ln_prior.min* and *ln_prior.max*.
    """
    outside = -np.inf

    for (parameter, value) in ln_prior.min.items():
        index = parameters_to_fit.index(parameter)
        if theta[index] < value:
            return outside

    for (parameter, value) in ln_prior.max.items():
        index = parameters_to_fit.index(parameter)
        if theta[index] > value:
            return outside

# Below we calculate prior probability based on x_caustic_in and x_caustic_out.
# This calculation assumes flat prior in (t_0, u_0, t_E, alpha), not in
# (x_caustic_in, x_caustic_out, t_caustic_in, t_caustic_out). If you want flat
# prior in the latter, then just replace following lines by "return 0".
    inside = event.model.parameters.uniform_caustic_sampling.jacobian(
        x_caustic_in=theta[parameters_to_fit.index('x_caustic_in')],
        x_caustic_out=theta[parameters_to_fit.index('x_caustic_out')])
    if inside == 0.:
        return outside
    else:
        return np.log(inside)


def ln_prob(
        theta, event, parameters_to_fit, print_models=False):
    """
    Log probability of the model - combines ln_prior() and ln_like().
    """
    ln_prior_ = ln_prior(theta, parameters_to_fit, event)
    if not np.isfinite(ln_prior_):
        return -np.inf

    ln_like_ = ln_like(theta, event, parameters_to_fit, print_models)
    if np.isnan(ln_like_):
        return -np.inf

    return ln_prior_ + ln_like_


def generate_random_parameters(parameters, starting, n, s=None, q=None):
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
    if 'x_caustic_in' in parameters and 'x_caustic_out' in parameters:
        sampling = mm.UniformCausticSampling(s=s, q=q)
        (x_in, x_out) = sampling.get_uniform_sampling(n)
        values[parameters.index('x_caustic_in')] = x_in
        values[parameters.index('x_caustic_out')] = x_out
    return np.array(values).T.tolist()


# Read config file.
config_file = "example_13.cfg"
config = configparser.ConfigParser()
config.optionxform = str  # So that "t_E" is not changed to "t_e".
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
k = {'comments': ['\\', '|']}
datasets = [mm.MulensData(file_name=f[0], phot_fmt=f[1], **k) for f in files]

# Generate starting values of parameters.
s = fixed_parameters.get('s', None)
q = fixed_parameters.get('q', None)
start = generate_random_parameters(
    parameters, starting, emcee_settings['n_walkers'], s=s, q=q)

# Setup Event instance that combines model and data.
par = dict(zip(parameters, start[0]))
par = {**par, **fixed_parameters}
my_model = mm.Model(par, coords=model_settings['coords'])
if 'methods' in model_settings:
    my_model.set_magnification_methods(model_settings['methods'])
if 'default_method' in model_settings:
    my_model.default_magnification_method = model_settings['default_method']
my_event = mm.Event(datasets=datasets, model=my_model)

# Prepare sampler.
n_dim = len(parameters)
print_models = other_settings.get('print_models', False)
args = (my_event, parameters, print_models)
sampler = emcee.EnsembleSampler(
    emcee_settings['n_walkers'], n_dim, ln_prob, args=args)

# Run sampler.
sampler.run_mcmc(start, emcee_settings['n_steps'])

# Parse results.
n_burn = emcee_settings['n_burn']
samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))
r_16 = np.percentile(samples, 16, axis=0)
r_50 = np.percentile(samples, 50, axis=0)
r_84 = np.percentile(samples, 84, axis=0)
print("\nFitted parameters:")
for i in range(n_dim):
    if parameters[i] == 'q':
        fmt = "{:} {:.7f} +{:.7f} -{:.7f}"
    else:
        fmt = "{:} {:.5f} +{:.5f} -{:.5f}"
    print(fmt.format(parameters[i], r_50[i], r_84[i]-r_50[i], r_50[i]-r_16[i]))
# We extract best model parameters and chi2 from the chain:
prob = sampler.lnprobability[:, n_burn:].reshape((-1))
best_index = np.argmax(prob)
best_chi2 = prob[best_index] / -0.5
best = samples[best_index, :]
print("\nSmallest chi2 model:")
print(*[repr(b) if isinstance(b, float) else b.value for b in best])
print(best_chi2)
for (best_, parameter) in zip(best, parameters):
    setattr(my_event.model.parameters, parameter, best_)

my_event.fit_fluxes()

# Expected results:
# t_0 ~ 2452848.06
# u_0 ~ 0.132
# t_E ~ 61.5
# alpha ~ 43.7
# It is possible that degenerate solution is found and then u_0 ~ -0.132 and
# alpha ~ 316.3
# You can inspect the output file and search for models similar to the above.
# The first one should appear within first 600 models calculated.
if 'x_caustic_in' in parameters:
    print(' t_0 = {:.5f}'.format(my_event.model.parameters.t_0))
    print(' u_0 = {:.5f}'.format(my_event.model.parameters.u_0))
    print(' t_E = {:.3f}'.format(my_event.model.parameters.t_E))
    print(' alpha = {:.2f}\n'.format(my_event.model.parameters.alpha))
print("chi2: ", my_event.get_chi2())  # Expected value: ~1655
