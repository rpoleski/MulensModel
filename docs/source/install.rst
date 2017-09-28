How to install?
===============

1. Make sure you have python with `astropy package`_ installed.
2. Download source code - either `recent release`_ or the current repository from `MulensModel github page`_ (green button on right).
3. Unpack the archive.
4. Add the path to the unpack directory to the ``PYTHONPATH``, e.g., in tcsh::

   setenv PYTHONPATH /home/USER_NAME/MulensModel-0.1.0/source\:$PYTHONPATH

in bash::

   export PYTHONPATH=/home/USER_NAME/MulensModel-0.1.0/source:$PYTHONPATH

5. Go to subdirecotry ``source/VBBL/`` and run ``make`` command. If it's not working and you're using Windows, then please run::

   gcc -lm -lstdc++ -fPIC -c VBBinaryLensingLibrary.cpp

   gcc -Wl,-soname,rapper -shared -o VBBinaryLensingLibrary_wrapper.so VBBinaryLensingLibrary_wrapper.cpp -lm -lstdc++ -fPIC VBBinaryLensingLibrary.o

6. Congratulations! You have MulensModel installed fully.

.. _astropy package: http://www.astropy.org/
.. _recent release: https://github.com/rpoleski/MulensModel/releases
.. _MulensModel github page: https://github.com/rpoleski/MulensModel
