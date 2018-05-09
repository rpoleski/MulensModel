"""
Fits PSPL model with parallax using EMCEE sampler.

"""
import os
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

from MulensModel import Event, Model, MulensData, MODULE_PATH


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


# Read the data
file_name = os.path.join(
    MODULE_PATH, "data/photometry_files", "starBLG234.6.I.218982.dat")
my_data = MulensData(file_name=file_name, add_2450000=True)

coords = "18:04:45.71 -26:59:15.2"

# Starting parameters:
params = dict()
params['t_0'] = 2453628.3
params['t_0_par'] = 2453628.
params['u_0'] = 0.37  # Change sign of u_0 to find the other solution.
params['t_E'] = 100.
params['pi_E_N'] = 0.
params['pi_E_E'] = 0.
my_model = Model(params, coords=coords)
my_event = Event(datasets=my_data, model=my_model)

# Which parameters we want to fit?
parameters_to_fit = ["t_0", "u_0", "t_E", "pi_E_N", "pi_E_E"]
# And remember to provide dispersions to draw starting set of points
sigmas = [0.01, 0.001, 0.1, 0.01, 0.01]

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
print("\nBest model:")
best = [my_event.best_chi2_parameters[p] for p in parameters_to_fit]
print(*[repr(b) if isinstance(b, float) else b.value for b in best])
print(my_event.best_chi2)

# Now let's plot 3 models
plt.figure()
model_0 = Model({'t_0': 2453628.29062, 'u_0': 0.37263, 't_E': 102.387105})
model_1 = Model(
    {'t_0': 2453630.35507, 'u_0': 0.488817, 't_E': 93.611301,
     'pi_E_N': 0.2719, 'pi_E_E': 0.1025, 't_0_par': params['t_0_par']},
    coords=coords)
model_2 = Model(
    {'t_0': 2453630.67778, 'u_0': -0.415677, 't_E': 110.120755,
     'pi_E_N': -0.2972, 'pi_E_E': 0.1103, 't_0_par': params['t_0_par']},
    coords=coords)
model_0.set_datasets([my_data])
model_1.set_datasets([my_data])
model_2.set_datasets([my_data])

t_1 = 2453200.
t_2 = 2453950.
plot_params = {'lw': 2.5, 'alpha': 0.3, 'subtract_2450000': True,
               't_start': t_1, 't_stop': t_2}

my_event.plot_data(subtract_2450000=True)
model_0.plot_lc(label='no pi_E', **plot_params)
model_1.plot_lc(label='pi_E, u_0>0', **plot_params)
model_2.plot_lc(label='pi_E, u_0<0', color='black', ls='dashed', **plot_params)

plt.xlim(t_1-2450000., t_2-2450000.)
plt.legend(loc='best')
plt.title('Data and 3 fitted models')
plt.show()
