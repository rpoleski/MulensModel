Basic Fitting Tutorial
======================

This tutorial shows how to fit basic point-source/point-lens model to 
the data. Similar material can be found in 
`example 2 <https://github.com/rpoleski/MulensModel/blob/master/examples/example_02_fitting.py>`_.

Importing data
--------------

First things first - we need to import some modules:

.. code-block:: python

   import os
   import matplotlib.pyplot as plt
   import MulensModel
   from MulensModel import Event, Fit, Model, MulensData, Utils

Then we import the data (downloaded together with the code) to 
the `MulensData class <https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html>`_:

.. code-block:: python

   file_name = os.path.join(MulensModel.MODULE_PATH, "data", 
                                                        "phot_ob08092_O4.dat")
   my_data = MulensData(file_name=file_name)
   print("{:} file was imported".format(file_name))

Plotting data
-------------

Next step would be plotting these data using matplotlib package:

.. code-block:: python

   plt.figure()
   plt.errorbar(my_data.time, my_data.mag, yerr=my_data.err_mag, fmt='.')
   plt.gca().invert_yaxis() # We need this to invert magnitude axis.
   plt.show()

From the plot we see that the peak was around JD' of 5380.0, 
the peak magnitude was about 0.8 mag brighter than baseline, 
and the event lasted dozens of days. 
We can turn these pieces of information into a very rough estimates of 
the event parameters:

.. code-block:: python

   t_0 = 5380.0
   u_0 = 0.5
   t_E = 20.0 # This is in days.

We guessed ``u_0 = 0.5`` based on the peak amplitude. The magnitude difference 
of 0.8 mag corresponds to flux ratio of slightly above 2. The magnification 
*A* and the impact parameter *u_0* are very approximately related via *A=1/u_0* 
so *u_0 = 0.5* should be a good choice. 

Preparing for fitting
---------------------

The rough estimates of the event parameters allow us to define 
a `Model <https://rpoleski.github.io/MulensModel/MulensModel.model.html>`_
and plot it 
together with data:

.. code-block:: python
   
   pspl_model = Model(t_0=t_0, u_0=u_0, t_E=t_E)
   pspl_model.set_datasets([my_data])
   pspl_model.plot_data()
   pspl_model.plot_lc()
   plt.show()

To associate a dataset with model we called Model.set_datasets() function. 
Its argument is a list of datasets. In the present case this list contains only 
a single dataset. If you have more datasets, then just include all of them
in the list passed to set_datasets(), e.g., 
``pspl_model.set_datasets([my_data, my_friends_data])``. 

The plot looks seems fine, i.e., the peak is more or less where it should be. 
Hence, we can use our rough estimates as a starting point for fitting 
procedure. 

You may want to learn more on plotting in MulensModel from 
`example 5 <https://github.com/rpoleski/MulensModel/blob/master/examples/example_05_MB08310.py>`_.

To fit the model parameters we will need to calculate chi^2. The chi^2 
calculation requires an instance of the 
`Event Class <https://rpoleski.github.io/MulensModel/MulensModel.event.htl>`_
that combines the Model and 
the MulensData:

.. code-block:: python
   
   my_event = Event(datasets=my_data, model=pspl_model)
   chi2_initial = my_event.get_chi2()
   print(my_event.model.parameters)
   print("give chi^2 of {:.2f}.".format(chi2_initial))

We have the ability to get the goodness of fit and it turn it into a function:

.. code-block:: python

   parameters_to_fit = ["t_0", "u_0", "t_E"]
   initial_guess = [t_0, u_0, t_E]

   def chi2_for_model(theta, event, parameters_to_fit):
       """for given event set attributes from parameters_to_fit (list of str) 
       to values from the theta list"""
       for (key, parameter) in enumerate(parameters_to_fit):
           setattr(event.model, parameter, theta[key])
       return event.get_chi2()

The chi2_for_model() function as a first argument has a sequence of 
float-type values. The second argument is an instance of the Event class. 
The third argument is a list that specifies the attributes of Event.model that 
will be changed. Note that the order of theta values and parameters_to_fit are 
the same. 

Fitting model parameters
------------------------

Ok, finally we can fit the parameters. Here we will use 
`the minimize() function <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html>`_ 
from SciPy subpackage optimize. We encourage you to 
try other fitting routines.

.. code-block:: python
   
   import scipy.optimize as op
   result = op.minimize(chi2_for_model, x0=initial_guess, 
           args=(my_event, parameters_to_fit), method='Nelder-Mead')

Fitting is done, so we can inspect the results. The function minimize() 
gives different output depending on method parameter. We will use just 
a few:

.. code-block:: python

   print("Fitting was successful? {:}".format(result.success))
   if not result.success:
       print(result.message)
   print("Function evaluations: {:}".format(result.nfev))
   print("The smallest function value: {:}".format(result.fun))
   print("for parameters: {:.5f} {:.4f} {:.3f}".format(*result.x.tolist()))

The best-fitting function parameters are stored in ``result.x``, which is 
of numpy.ndarray type. To have a nice output, we converted them to a list. 
Let's plot two different models:

.. code-block:: python

   # Initial model:
   pspl_model.t_0 = t_0
   pspl_model.u_0 = u_0
   pspl_model.t_E = t_E
   pspl_model.plot_lc(label='initial', c='red')
   # Best fitting model:
   pspl_model.t_0 = result.x[0]
   pspl_model.u_0 = result.x[1]
   pspl_model.t_E = result.x[2]
   pspl_model.plot_lc(label='fitted')
   # Finally: data, legend, and show the plot:
   pspl_model.plot_data()
   plt.legend(loc='best')
   plt.show()

If you zoom-in on the peak, you will easily see that the fitted model is 
much better. 

Congratulations! You have fitted the model to the data.

Excercise
---------

Try using different optimization routine, starting point, 
or apply constraints on the fit. If 
`the minimize() function <https://docs.scipy.org/doc/scipy/reference/optimize.html>`_ 
is now your favourite fitting routine, then still you can call it differently. 
Try changing ``method`` parameter to one of: 
'Powell', 'CG', 'BFGS', 'TNC', 'COBYLA'.

