.. MulensModel documentation master file, created by
   sphinx-quickstart on Wed Aug 23 15:22:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MulensModel's documentation!
=======================================

MulensModel is package for modeling microlensing (or :math:`\\mu`-lensing) 
events. MulensModel can generate a microlensing light curve for a given set of 
microlensing parameters, fit that light curve to some data, and return a chi2 
value. That chi2 can then be input into an arbitrary likelihood function to 
find the best fit parameters.

Main MulensModel classes are `MulensData`_, `Model`_, `ModelParameters`_, 
and `Event`_. For descriptions of all classes see :ref:`modindex`.

.. _MulensData: https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html
.. _Model: https://rpoleski.github.io/MulensModel/MulensModel.model.html
.. _ModelParameters: https://rpoleski.github.io/MulensModel/MulensModel.modelparameters.html
.. _Event: https://rpoleski.github.io/MulensModel/MulensModel.event.html

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Installation <install>
   Tutorial <tutorial>
   Fitting tutorial <tutorial_fit_pspl>
   Parallax fitting tutorial <tutorial_fit_pi_E>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

