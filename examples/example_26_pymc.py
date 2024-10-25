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

import torch as pt
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
parameters_to_fit = ['t_0']


def ln_like(t_0):
    """ likelihood function """
    print(type(t_0))
    event.model.parameters.t_0 = t_0
    #event.model.parameters.u_0 = u_0
    #event.model.parameters.t_E = t_E
    #event.model.parameters.t_star = t_star

    chi2 = event.get_chi2()

    return -0.5 * chi2


class LogLike(pt.Op):

    def make_node(self, t_0) -> Apply:
        # Convert inputs to tensor variables
        t_0 = pt.as_tensor(t_0)
        inputs = [t_0]
        # Define output type, in our case a vector of likelihoods
        # with the same dimensions and same data type as data
        # If data must always be a vector, we could have hard-coded
        # outputs = [pt.vector()]
        outputs = float

        # Apply is an object that combines inputs, outputs and an Op (self)
        return Apply(self, inputs, outputs)

    def perform(self, node: Apply, inputs: list[np.ndarray], outputs: list[list[None]]) -> None:
        # This is the method that compute numerical output
        # given numerical inputs. Everything here is numpy arrays
        t_0 = inputs  # this will contain my variables

        # call our numpy log-likelihood function
        loglike_eval = ln_like(t_0)

        # Save the result in the outputs list provided by PyTensor
        # There is one list per output, each containing another list
        # pre-populated with a `None` where the result should be saved.
        outputs[0][0] = np.asarray(loglike_eval)


basic_model = pm.Model()
with basic_model:
    # Priors for unknown model parameters
    t_0 = pm.Normal('t_0', mu=2454656.4, sigma=0.001)
    #u_0 = pm.Normal('u_0', mu=0.003, sigma=0.0003)
    #t_E = pm.Normal('t_E', mu=11.14, sigma=0.05)
    #t_star = pm.Normal('t_star', mu=0.055, sigma=0.001)

    pm.CustomDist('likelihood', t_0, u_0, t_E, t_star, logp=ln_like)

with basic_model:
    # draw 1000 posterior samples
    idata = pm.sample()

az.summary(idata)
az.plot_trace(idata, combine=True)
