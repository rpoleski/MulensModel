How to install?
===============

1. Make sure you have python with `astropy package`_ installed. Note that newest astropy release requires python 3.X, so if you only have python 2.X, than install older version of astropy, e.g., 1.3.3. 
2. Download source code - either `recent release`_ or the current repository from `MulensModel github page`_ (green button on right).
3. Unpack the archive.
4. Add the path to the unpack directory to the ``PYTHONPATH``, e.g., if you extracted archive in your home directory (``/home/USER_NAME/``) in tcsh::

   setenv PYTHONPATH /home/USER_NAME/MulensModel-0.2.1/source\:$PYTHONPATH

in bash::

   export PYTHONPATH=/home/USER_NAME/MulensModel-0.2.1/source:$PYTHONPATH

In order to have this command invoked every time you open the terminal, please add this command to ``~/.cshrc`` or ``~/.bashrc`` file.
5. Go to subdirecotry ``source/VBBL/`` and run ``make`` command. If it's not working and you're using Windows, then please run::

   gcc -lm -lstdc++ -fPIC -c VBBinaryLensingLibrary.cpp

   gcc -Wl,-soname,rapper -shared -o VBBinaryLensingLibrary_wrapper.so VBBinaryLensingLibrary_wrapper.cpp -lm -lstdc++ -fPIC VBBinaryLensingLibrary.o

6. Congratulations! You have MulensModel installed fully.

.. _astropy package: http://www.astropy.org/
.. _recent release: https://github.com/rpoleski/MulensModel/releases
.. _MulensModel github page: https://github.com/rpoleski/MulensModel
