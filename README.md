# MulensModel

MulensModel is package for modelling microlensing (or $\mu$-lensing) events. 

It is still under development. 

MulensModel can generate a microlensing light curve for a given set of microlensing parameters, fit that light curve to some data, and return a chi2 value. That chi2 can then be input into an arbitrary likelihood function to find the best fit parameters.

There is no tutorial yet, but we have a few useful resources:

* [Examples on how to use the code](examples/):
  * [Example 01](examples/example_01_models.py) -- plot simple point-source/point-lens (PSPL) model and model with planetary lens,
  * [Example 02](examples/example_02_fitting.py) -- fit PSPL model to the data using scipy.optimize.minimize(),
  * [Example 03](examples/example_03_mulenssystem.py) -- define PSPL model using physical properties and plot the resulting magnification curve,
  * [Example 04](examples/example_04_einsteinring.py) -- calculate the Einstein ring size for a grid of lens masses and distances,
  * [Example 05](examples/example_05_MB08310.py) -- plot multiple datasets for a single model, plot the residuals, and do this both in magnitude and magnification spaces,
* [Instructions on getting satellite positions](documents/Horizons_manual.md)

More will be added soon.

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

---

file revised Apr 2017

