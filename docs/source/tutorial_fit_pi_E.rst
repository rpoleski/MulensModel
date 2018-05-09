Microlensing Parallax Fitting Tutorial
======================================

Here you will learn how to fit the observed light-curve with a model that 
has a point lens, point source, but their relative motion is not rectilinear 
and includes annual microlensing parallax effect. If you haven't yet looked 
at the other tutorials, then you may want to start there and come back here 
later.

We need some data to fit. Let's look at some long-timescale event from 
`Wyrzykowski et al. (2015) 
<http://adsabs.harvard.edu/abs/2015ApJS..216...12W>`_. As an example 
we can take event with ID 3291 (star BLG234.6.I.218982; also named 
OGLE-2005-BLG-086 on `OGLE EWS website
<http://ogle.astrouw.edu.pl/ogle4/ews/ews.html>`_). We can look at the plot of
model without parallax: `here 
<ftp://ftp.astrouw.edu.pl/ogle/ogle3/blg_tau/PLOTS/starBLG234.6.I.218982.dat.png>`_.
We see that model does not fit data very well - there are trends in 
the residuals around 3400 and 3800. 

Imports and settings
--------------------

We start by importing python modules:

.. code-block:: python

   import os
   import numpy as np
   import emcee
   import matplotlib.pyplot as plt
   import MulensModel
   from MulensModel import Event, Model, MulensData

We're using EMCEE package for fitting. You can download it from
`<http://dfm.io/emcee/current/>`_ in case you don't have it yet. Then we
import the data and set the event coordinates:

.. code-block:: python

   file_name = os.path.join(MulensModel.MODULE_PATH, "data", 
       "photometry_files", "starBLG234.6.I.218982.dat")
   my_data = MulensData(file_name=file_name, add_2450000=True)
   coords = "18:04:45.71 -26:59:15.2"

Note that *add_2450000=True* is very important. The file has time vector 
of HJD-2450000 for convenience and it's fine as long as we don't fit 
the annual parallax. He have to set *MulensData* time vector to
full HJD because this is used to calculate the position of Earth 
relative to Sun and we don't want this to be off by around 6700 years. 

We set the starting values of the parameters:

.. code-block:: python

   params = dict()
   params['t_0'] = 2453628.3
   params['t_0_par'] = 2453628.
   params['u_0'] = 0.37
   params['t_E'] = 100.
   params['pi_E_N'] = 0.
   params['pi_E_E'] = 0.
   my_model = Model(params, coords=coords)
   my_event = Event(datasets=my_data, model=my_model)

We set the parameter reference time (*t_0_par*) for rounded value of *t_0*.
This is common approach. If you don't set *t_0_par*, then fitting will be 
slower, because Earth positions will be re-calculated for every model. 

Further we need to specifies which parameters we want to fit and also 
specify dispersions in starting points. We choose the latter to be some 
small values:

.. code-block:: python

   parameters_to_fit = ["t_0", "u_0", "t_E", "pi_E_N", "pi_E_E"]
   sigmas = [0.01, 0.001, 0.1, 0.01, 0.01]

Some more EMCEE settings - number of walkers, steps, and burn-in steps. Also
the list of starting points:

.. code-block:: python

   n_dim = len(parameters_to_fit)
   n_walkers = 40
   n_steps = 500
   n_burn = 150
   start_1 = [params[p] for p in parameters_to_fit]
   start = [start_1 + np.random.randn(n_dim) *  sigmas
           for i in range(n_walkers)]

We need one more important piece of information - the function that 
computes the logarithm of (unnormalized) probability. We split it into
three separate functions for clarity:

.. code-block:: python

   def ln_like(theta, event, parameters_to_fit):
       """ likelihood function """
       for key, val in enumerate(parameters_to_fit):
           setattr(event.model.parameters, val, theta[key])
       return -0.5 * event.get_chi2()

.. code-block:: python
   
   def ln_prior(theta, parameters_to_fit):
       """priors - we only reject obviously wrong models"""
       if theta[parameters_to_fit.index("t_E")] < 0.:
           return -np.inf
       return 0.0

.. code-block:: python

   def ln_prob(theta, event, parameters_to_fit):
       """ combines likelihood and priors"""
       ln_prior_ = ln_prior(theta, parameters_to_fit)
       if not np.isfinite(ln_prior_):
           return -np.inf
       ln_like_ = ln_like(theta, event, parameters_to_fit)
       if np.isnan(ln_like_): 
           return -np.inf
       return ln_prior_ + ln_like_
   
Running the sampler
-------------------

Ok, we're ready to run EMCEE:

.. code-block:: python

   sampler = emcee.EnsembleSampler(
       n_walkers, n_dim, ln_prob, args=(my_event, parameters_to_fit))
   sampler.run_mcmc(start, n_steps)
   samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))

And now we're ready to look at the results and best-fitted model:

.. code-block:: python

   results = np.percentile(samples, [16, 50, 84], axis=0)
   print("Fitted parameters:")
   form = "{:.5f} {:.5f} {:.5f}"
   for i in range(n_dim):
       r = results[1, i]
       print(form.format(r, results[2, i]-r, r-results[0, i]))
   print("\nBest model:")    
   best = [my_event.best_chi2_parameters[p] for p in parameters_to_fit]
   print(*[repr(b) if isinstance(b, float) else b.value for b in best])
   print(my_event.best_chi2)

I hope you got (u_0, t_E, pi_E_N, pi_E_E) of around
(0.44, 95, 0.21, 0.10) and chi^2 of 949.5. 

At this point you may want to say that the fit is done at this point.
But it's not! We have to check for degenerate solution. We're fitting single
lens model, hence, the search for degenerate solution is easy and it's enough
to start with negative u_0. 

Now you have time to do the second fit...

Ok, I hope you got (u_0, t_E, pi_E_N, pi_E_E) of
(-0.41, 110, -0.30, 0.11) and chi^2 of 947.0. The difference between
the two solutions is small in chi^2 - they are degenerate. And u_0<0 fits
data slightly better. It turned out that the second fit was very important!

Plotting
--------

Let's make a nice plot! 

I provide model parameters below. Here is how it goes:

.. code-block:: python

   plt.figure()
   model_0 = Model({'t_0': 2453628.29062, 'u_0': 0.37263,
           't_E': 102.387105})
   model_1 = Model({'t_0': 2453630.35507, 'u_0': 0.488817,
           't_E': 93.611301, 'pi_E_N': 0.2719, 'pi_E_E': 0.1025,
           't_0_par': params['t_0_par']}, coords=coords)
   model_2 = Model({'t_0': 2453630.67778, 'u_0': -0.415677,
           't_E': 110.120755, 'pi_E_N': -0.2972, 'pi_E_E': 0.1103,
           't_0_par': params['t_0_par']}, coords=coords)
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
   model_2.plot_lc(label='pi_E, u_0<0', color='black', ls='dashed',
           **plot_params)
   
   plt.xlim(t_1-2450000., t_2-2450000.)
   plt.legend(loc='best')
   plt.title('Data and 3 fitted models')
   plt.show()

I hope you see that parallax models are better than the non-parallax model.
If not, then zoom-in around epoch 3800. The non-parallax model has chi^2
higher by about 400.

Slightly modified source code from this tutorial is
`example 6 
<https://github.com/rpoleski/MulensModel/blob/master/examples/example_06_fit_parallax_EMCEE.py>`_.
Additionally, `example 7 
<https://github.com/rpoleski/MulensModel/blob/master/examples/example_07_fit_parallax_MN.py>`_ 
shows how to fit parallax model using MultiNest instead of EMCEE algorithm.  
Note that a single run of MultiNest finds two degenerate modes and reports 
properties of both of them.  


Exercise
--------

As an exercise you may try to fit other events from 
`Wyrzykowski et al. (2015)`_. It's best to start with long events, that have bright sources, and small impact parameters.

