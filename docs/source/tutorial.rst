Very Basic Tutorial
===================

The main MulensModel features are in classes `Model`_, `MulensData`_,
and `Event`_.

This is a very simple tutorial showing how you might use those classes. It is also available (somewhat expanded) as a Jupyter `notebook`_.

.. _Model: https://rpoleski.github.io/MulensModel/MulensModel.model.html
.. _Event: https://rpoleski.github.io/MulensModel/MulensModel.event.html
.. _MulensData: https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html
.._notebook: https://www.github.com/rpoleski/MulensModel/examples/MulensModelTutorial.ipynb

Defining a Model
----------------

The most basic thing to do is to define a microlensing model. For example, you could define a point lens model as follows:

``my_pspl_model = MulensModel.Model(t_0=2452848.06, u_0=0.133, t_E=61.5)``

Or a model with 2-bodies:

``my_1S2L_model = MulensModel.Model(t_0=2452848.06, u_0=0.133, t_E=61.5, rho=0.00096, q=0.0039, s=1.120, alpha=43.8)``

(by default alpha is in degrees, but you could explicitly specify radians)

Then, you might plot those models:

``import matplotlib.pyplot as pl``

``pl.figure()``

``my_pspl_model.plot_magnification(t_range=[2452810,2452890], subtract_2450000=True, color='red', linestyle=':')``

``my_1S2L_model.plot_magnification(t_range=[2452810,2452890], subtract_2450000=True, color='black')``

``pl.show()``

Introducing Data
----------------

Suppose you also had some data you want to import:

``OGLE_data = MulensModel.MulensData('../data/OB03235/OB03235_OGLE.tbl.txt', commenst=['\\','|'])``

``MOA_data = MulensModel.MulensData('../data/OB03235/OB03235_MOA.tbl.txt', phot_fmt='flux', comments=['\\','|')``

Combining Data with a Model
---------------------------

Now suppose you wanted to combine the two together:

``my_event = MulensModel.Event(datasets=[OGLE_data, MOA_data], model=my_1S2L_model)``

And you wanted to plot the result:

``my_event.plot_model(t_range=[2452810, 2452890], subtract_2450000=True, color='black')``

``my_event.plot_data(subtract_2450000=True)``

``pl.show()``

This fits for the fluxes so that the model and data are all on the
flux scale set by the first dataset. It does NOT fit for the best
microlensing parameters. If you wanted to know how good the fit is, you can get the chi2:

``print(my_event.get_chi2())``

If you want to optimize that chi2, we leave it up to you to determine the best method for doing this.