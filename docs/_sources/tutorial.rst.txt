Very Basic Tutorial
===================

The main MulensModel features are in classes `Model`_, `MulensData`_,
and `Event`_.

This is a very simple tutorial showing how you might use those classes. 
It is also available (somewhat expanded) as a Jupyter `notebook`_. 
Please note that MulensModel is written in Python3. If you're using Python2.X, 
then start by adding ``from __future__ import print_function`` at the begin 
of your codes and be advised that we don't guarantee that everything will work. 

.. _Model: https://rpoleski.github.io/MulensModel/MulensModel.model.html
.. _Event: https://rpoleski.github.io/MulensModel/MulensModel.event.html
.. _MulensData: https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html
.. _notebook: https://github.com/rpoleski/MulensModel/blob/master/examples/MulensModelTutorial.ipynb   

This example shows OGLE-2003-BLG-235/MOA-2003-BLG-53, the first
microlensing planet. See `Bond et al. (2004) 
<https://ui.adsabs.harvard.edu/abs/2004ApJ...606L.155B/abstract>`_.
The data were downloaded from the `NASA Exoplanet Archive
<https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname=OGLE-2003-BLG-235L+b&type=CONFIRMED_PLANET>`_.

Defining a Model
----------------

The most basic thing to do is to define a microlensing model. For example, you could define a point lens model as follows:

.. code-block:: python

   import MulensModel as mm
   my_pspl_model = mm.Model({'t_0': 2452848.06, 'u_0': 0.133, 't_E': 61.5})

Or a model with 2-bodies:

.. code-block:: python
   
   my_1S2L_model = mm.Model({'t_0': 2452848.06, 'u_0': 0.133, 
        't_E': 61.5, 'rho': 0.00096, 'q': 0.0039, 's': 1.120, 
        'alpha': 223.8})

(by default alpha is in degrees, but you could explicitly specify radians)

Since rho is set, define a time range and method for finite source 
effects:

.. code-block:: python

   my_1S2L_model.set_magnification_methods([2452833., 'VBBL', 2452845.])

Then, you might plot those models:

.. code-block:: python
   
   import matplotlib.pyplot as plt
   my_pspl_model.plot_magnification(t_range=[2452810, 2452890], 
       subtract_2450000=True, color='red', linestyle=':')
   my_1S2L_model.plot_magnification(t_range=[2452810, 2452890], 
       subtract_2450000=True, color='black')
   plt.show()

Introducing Data
----------------

Suppose you also had some data you want to import:

.. code-block:: python

   import os
   path = os.path.join(mm.DATA_PATH, 'photometry_files', 'OB03235')
   OGLE_data = mm.MulensData(
        file_name=os.path.join(path, 'OB03235_OGLE.tbl.txt'),
        comments=['\\', '|'])
   MOA_data = mm.MulensData(
        file_name=os.path.join(path, 'OB03235_MOA.tbl.txt'),
        comments=['\\', '|'], phot_fmt='flux')

Combining Data with a Model
---------------------------

Now suppose you wanted to combine the two together:

.. code-block:: python

   my_event = mm.Event(datasets=[OGLE_data, MOA_data], 
       model=my_1S2L_model)

And you wanted to plot the result:

.. code-block:: python
   
   my_event.plot_model(t_range=[2452810, 2452890], subtract_2450000=True, 
       color='black')
   my_event.plot_data(subtract_2450000=True)
   plt.xlim(2810, 2890)
   plt.ylim(19.25, 16.6)
   plt.show()

This fits for the fluxes so that the model and data are all on the
flux scale set by the first dataset. It does NOT fit for the best
microlensing parameters. If you wanted to know how good the fit is, you can get the chi2:

.. code-block:: python
   
   print(my_event.get_chi2())

If you want to optimize that chi2, we leave it up to you to determine the best method for doing this.

