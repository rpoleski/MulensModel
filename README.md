# MulensModel

<dl>MulensModel is package for modeling microlensing (or &mu;-lensing) 
events. </dl>

It is still under development. [Latest release: 0.2.1](https://github.com/rpoleski/MulensModel/releases/tag/0.2.1)

MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best fit parameters.

A few useful resources:

* [Basic usage tutorial](https://rpoleski.github.io/MulensModel/tutorial.html),
* [Fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html),
* [Examples on how to use the code](examples/):
  * [Example 01](examples/example_01_models.py) -- plot simple point-source/point-lens (PSPL) model and model with planetary lens,
  * [Example 02](examples/example_02_fitting.py) -- fit PSPL model to the data using scipy.optimize.minimize(),
  * [Example 03](examples/example_03_mulenssystem.py) -- define PSPL model using physical properties and plot the resulting magnification curve,
  * [Example 04](examples/example_04_einsteinring.py) -- calculate the Einstein ring size for a grid of lens masses and distances,
  * [Example 05](examples/example_05_MB08310.py) -- plot multiple datasets for a single model, plot the residuals, and do this both in magnitude and magnification spaces,
* [Instructions on getting satellite positions](documents/Horizons_manual.md)

More will be added soon.

[Documentation](https://rpoleski.github.io/MulensModel/)

If you want to learn more about microlensing, please visit [Microlensing Source website](http://microlensing-source.org/).

Currently, MulensModel supports:
* Lens Systems: Point Lens, Binary Lens
* Source Stars: Single source
* Effects: finite source (1-parameter)

Under Development:
* Effects: finite source (2-parameter), parallax (satellite, annual)

Future Development:
* Source Stars: Binary source, xallarap
* Effects: parallax (topographic)

### How to install?

1. Make sure you have python with [astropy package](http://www.astropy.org/) installed.
2. Download source code - either [recent release](https://github.com/rpoleski/MulensModel/releases) or current repository using green button above. 
3. Unpack the archive.
4. Add the path to the unpack directory to the PYTHONPATH, e.g., if you've extracted the archive in your home directory (``/home/USER_NAME/``) in tcsh:
```
setenv PYTHONPATH /home/USER_NAME/MulensModel-0.2.1/source\:$PYTHONPATH
```
in bash:
```
export PYTHONPATH=/home/USER_NAME/MulensModel-0.2.1/source:$PYTHONPATH
```
In order to have this command invoked every time you open the terminal, please add this command to ``~/.cshrc`` or ``~/.bashrc`` file.

5. Go to subdirecotry ```source/VBBL/``` and run ```make``` command. If it's not working and you're using Windows, then please run:
```
gcc -lm -lstdc++ -fPIC -c VBBinaryLensingLibrary.cpp
gcc -Wl,-soname,rapper -shared -o VBBinaryLensingLibrary_wrapper.so VBBinaryLensingLibrary_wrapper.cpp -lm -lstdc++ -fPIC VBBinaryLensingLibrary.o
```
6. Congratulations! You have MulensModel installed fully. 

---
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

file revised Nov 2017

