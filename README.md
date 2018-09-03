# MulensModel

<dl>MulensModel is package for modeling microlensing (or &mu;-lensing) 
events. </dl>

[Latest release: 1.4.0](https://github.com/rpoleski/MulensModel/releases/latest) and we're working on further developing the code.

MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best fit parameters.

Are you interested in [Microlensing Data Analysis Challenge](http://microlensing-source.org/data-challenge/)? For more details on using MulensModel in data challenge, see [this document](documents/data_challenge.md).

A few useful resources:

* [Basic usage tutorial](https://rpoleski.github.io/MulensModel/tutorial.html),
* [Fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html),
* [Microlensing parallax fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pi_E.html),
* [Examples on how to use the code](examples/):
  * [Example 01](examples/example_01_models.py) -- plot simple point-source/point-lens (PSPL) model and model with planetary lens,
  * [Example 02](examples/example_02_fitting.py) -- fit PSPL model to the data using scipy.optimize.minimize(),
  * [Example 03](examples/example_03_mulenssystem.py) -- define PSPL model using physical properties and plot the resulting magnification curve,
  * [Example 04](examples/example_04_einsteinring.py) -- calculate the Einstein ring size for a grid of lens masses and distances,
  * [Example 05](examples/example_05_MB08310.py) -- plot multiple datasets for a single model, plot the residuals, and do this both in magnitude and magnification spaces,
  * [Example 06](examples/example_06_fit_parallax_EMCEE.py) -- fit parallax model using EMCEE,
  * [Example 07](examples/example_07_fit_parallax_MN.py) -- fit parallax model using MultiNest,
  * [Example 08](examples/example_08_planet_grid_fitting.ipynb) -- shows how to fit simulated WFIRST light curve with planetary model,
  * [Example 09](examples/example_09_gradient_fitting.py) -- fit point lens model using chi^2 gradient,
  * [Example 10](examples/example_10_fitting_and_fluxes.py) -- fit model and extract posterior fluxes.
* [Instructions on getting satellite positions](documents/Horizons_manual.md)

[Documentation](https://rpoleski.github.io/MulensModel/) includes description of input and output of every function. 

If you want to learn more about microlensing, please visit [Microlensing Source website](http://microlensing-source.org/).

Currently, MulensModel supports:
* Lens Systems: point lens or binary lens.
* Source Stars: single source or binary source.
* Effects: finite source (1-parameter), parallax (satellite & annual), binary lens orbital motion, different parametrizations of microlensing models.

Need more? Open [an issue](https://github.com/rpoleski/MulensModel/issues) or send us an e-mail. 

Are you using MulensModel for scientific research? Please give us credit by citing the [paper](http://adsabs.harvard.edu/abs/2018arXiv180301003P) and [ASCL reference](http://ascl.net/1803.006).

## How to install?

Download the source code and run:
```
python setup.py install
```
MulensModel requires some standard packages plus [astropy package](http://www.astropy.org/). To make sure that you have everything that's needed, just run:
```
pip install -r requirements.txt
```
If you have any problems, please contact the authors and we will try to help.

---
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

file revised Sep 2018

