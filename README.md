# MulensModel

<dl>MulensModel is package for modeling microlensing (or &mu;-lensing) 
events. </dl>

[Latest release: 1.14.0](https://github.com/rpoleski/MulensModel/releases/latest) and we're working on further developing the code.

MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best fit parameters.

A few useful resources:

* [Basic usage tutorial](https://rpoleski.github.io/MulensModel/tutorial.html),
* [Fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pspl.html),
* [Microlensing parallax fitting tutorial](https://rpoleski.github.io/MulensModel/tutorial_fit_pi_E.html),
* [**Explanation of microlensing parameters**](documents/parameter_names.pdf)
* [**Explanation of methods used for calculating magnification**](documents/magnification_methods.pdf)
* [Examples on how to use the code](examples/):
  * [Example 01](examples/example_01_models.py) - plot simple point-source/point-lens (PSPL) model and model with planetary lens,
  * [Example 02](examples/example_02_fitting.py) - fit PSPL model to the data using scipy.optimize.minimize(),
  * [Example 03](examples/example_03_mulenssystem.py) - define PSPL model using physical properties and plot the resulting magnification curve,
  * [Example 04](examples/example_04_einsteinring.py) - calculate the Einstein ring size for a grid of lens masses and distances,
  * [Example 05](examples/example_05_MB08310.py) - plot multiple datasets for a single model, plot the residuals, and do this both in magnitude and magnification spaces,
  * [Example 06](examples/example_06_fit_parallax_EMCEE.py) - fit parallax model using EMCEE,
  * [Example 07](examples/example_07_fit_parallax_MN.py) - fit parallax model using MultiNest,
  * [Example 08](examples/example_08_planet_grid_fitting.ipynb) - shows how to fit simulated WFIRST light curve with planetary model,
  * [Example 09](examples/example_09_gradient_fitting.py) - fit point lens model using chi^2 gradient,
  * [Example 10](examples/example_10_fitting_and_fluxes.py) - fit model and extract posterior fluxes, **use [config file](examples/example_10.cfg) to pass all parameters**,
  * [Example 11](examples/example_11_binary_source.py) - simulate and fit binary source event,
  * [Example 12](examples/example_12_fit_satellite_parallax_EMCEE.py) - fit parallax model to ground and satellite data, plot models and trajectories at the end,
  * [Example 13](examples/example_13_caustic_sampling.py) - fit planetary event using caustic entrance and exit epochs as parameters (uses [config file](examples/example_13.cfg)),
  * [Example 14](examples/example_14_caustic_plotting.py) - plot caustic using standard method and uniform sampling,
  * [Example 15](examples/example_15_fitting.py) - fitting binary lens model with many options - use config file for [ob05390](examples/example_15_ob05390_v1.cfg) or [mb07192](examples/example_15_mb07192_v1.cfg); settings are read by [this file](examples/example_15_read.py).
* [Instructions on getting satellite positions](documents/Horizons_manual.md).

[Documentation](https://rpoleski.github.io/MulensModel/) includes description of input and output of every function. 

If you want to learn more about microlensing, please visit [Microlensing Source website](http://microlensing-source.org/).

Currently, MulensModel supports:
* Lens Systems: point lens or binary lens.
* Source Stars: single source or binary source.
* Effects: finite source (1-parameter), parallax (satellite & annual), binary lens orbital motion, different parametrizations of microlensing models.

Need more? Open [an issue](https://github.com/rpoleski/MulensModel/issues) or send us an e-mail. 

Are you using MulensModel for scientific research? Please give us credit by citing the [paper published in "Astronomy and Computing"](https://ui.adsabs.harvard.edu/abs/2019A%26C....26...35P/abstract) and [ASCL reference](http://ascl.net/1803.006). For arXiv see [link](https://arxiv.org/abs/1803.01003).

## How to install?

Download the source code and run:
```
python setup.py install
```
MulensModel requires some standard packages plus [astropy package](http://www.astropy.org/). To make sure that you have everything that's needed, just run:
```
pip install -r requirements.txt
```
Alternatively, you can run makefiles: go to `source/VBBL/` and run `make`, then go to `source/AdaptiveContouring/` and do the same. Then and add the path `MulensModel/source` to your `PYTHONPATH`. If you have any problems, please contact the authors and we will try to help.

If you have **problems with installing or running MulensModel on MacOS**, please see notes [here](documents/macos_install.md).

---
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

file revised Oct 2019

