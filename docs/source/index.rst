.. MulensModel documentation master file, created by
   sphinx-quickstart on Wed Aug 23 15:22:44 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MulensModel's documentation!
=======================================

MulensModel is package for modeling microlensing (or :math:`\mu`-lensing) 
events. MulensModel can generate a microlensing light curve for a given set of 
microlensing parameters, fit that light curve to some data, and return a chi2 
value. That chi2 can then be input into an arbitrary likelihood function to 
find the best fit parameters.

Main MulensModel features are in classes `Model`_, `Event`_, and 
`MulensData`_. 

How to install?
===============

1. Make sure you have python with `astropy package`_ installed.
2. Download source code - either `recent release`_ or the current repository 
from `MulensModel github page`_ (green button on right).
3. Unpack the archive.
4. Add the path to the unpack directory to the *PYTHONPATH*, e.g., in tcsh::

   setenv PYTHONPATH /home/USER_NAME/MulensModel-0.1.0/source\:$PYTHONPATH

in bash::

   export PYTHONPATH=/home/USER_NAME/MulensModel-0.1.0/source:$PYTHONPATH

5. Go to subdirecotry *source/VBBL/* and run *make* command. If it's not working and you're using Windows, then please run::

   gcc -lm -lstdc++ -fPIC -c VBBinaryLensingLibrary.cpp

   gcc -Wl,-soname,rapper -shared -o VBBinaryLensingLibrary_wrapper.so VBBinaryLensingLibrary_wrapper.cpp -lm -lstdc++ -fPIC VBBinaryLensingLibrary.o

6. Congratulations! You have MulensModel installed fully.

.. _astropy package: http://www.astropy.org/
.. _recent release: https://github.com/rpoleski/MulensModel/releases
.. _MulensModel github page: https://github.com/rpoleski/MulensModel
.. _Model: https://rpoleski.github.io/MulensModel/MulensModel.model.html
.. _Event: https://rpoleski.github.io/MulensModel/MulensModel.event.html
.. _MulensData: https://rpoleski.github.io/MulensModel/MulensModel.mulensdata.html

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

