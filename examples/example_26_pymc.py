"""
Fit a PSPL model to data using PyMC. Based on Example 10.
"""
import sys
try:
    import pymc as pm
except ImportError as err:
    print(err)
    print("\npymc could not be imported.")
    print("Get it and re-run the script")
    sys.exit(1)

import arviz as az
import numpy as np
import MulensModel as mm

# JCY: For now, we're going to ignore bad data
file_names = ['../data/photometry_files/MB08310/MOA_0300089_PLC_007.tbl',
              '../data/photometry_files/MB08310/CTIO_H_0300089_PLC_004.tbl',
              '../data/photometry_files/MB08310/Bron_0300089_PLC_002.tbl']
kwargs = {'comments': ["\\", "|"]}
data = [mm.MulensData(file_name=name, **kwargs) for name in file_names]

parameters = {'t_0': 2454656.4, 'u_0': 0.003, 't_E': 11.14, 't_star': 0.055}
model = mm.Model(parameters)
model.set_magnification_methods(
    [2454656.25, 'finite_source_uniform_Gould94', 2454656.55])

event = mm.Event(datasets=data, model=model)

# New pymc stuff:
parameters_to_fit = ['t_0', 'u_0', 't_E', 't_star']
basic_model = pm.Model()

observed_fluxes = None
observed_errors = None
for dataset in event.datasets:
    if observed_fluxes is None:
        observed_fluxes = dataset.flux
        observed_errors = dataset.err_flux
    else:
        observed_fluxes = np.hstack([observed_fluxes, dataset.flux])
        observed_errors = np.hstack([observed_errors, dataset.err_flux])

with basic_model:
    # Priors for unknown model parameters
    t_0 = pm.Normal('t_0', mu=2454656.4, sigma=0.001)
    u_0 = pm.Normal('u_0', mu=0.003, sigma=0.0003)
    t_E = pm.Normal('t_E', mu=11.14, sigma=0.05)
    t_star = pm.Normal('t_star', mu=0.055, sigma=0.001)

    event.model.parameters.t_0 = t_0
    event.model.parameters.u_0 = u_0
    event.model.parameters.t_E = t_E
    event.model.parameters.t_star = t_star
    event.fit_fluxes()

    # Expected value of outcome
    mod_fluxes = None
    for fit in event.fits:
        if mod_fluxes is None:
            mod_fluxes = fit.get_model_fluxes()
        else:
            mod_fluxes = np.hstack([mod_fluxes, fit.get_model_fluxes()])

    # Likelihood (sampling distribution) of observations
    Y_obs = pm.Normal("Y_obs", mu=mod_fluxes, sigma=observed_errors, observed=observed_fluxes)

with basic_model:
    # draw 1000 posterior samples
    idata = pm.sample()

az.summary(idata)
az.plot_trace(idata, combine=True)
