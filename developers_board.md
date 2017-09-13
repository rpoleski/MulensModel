## current goals:
1. code review, 
2. instantaneous orbital motion of binary lens - use case and unit test,
3. docstrings for Spinx in all public methods of Event, Model, and MulensData classes,
4. reasonable website with documentation


## Specific tasks to be performed
(__boldfaced__ correspond to this month goals; _italics_ mark task useful for data challenge; try to put important stuff at the top)

* PSPL manual
* __correct documentation for existing code - start with Model, MulensData, and Event__
* anything from use cases that does not work yet -- see TODO.md file
* __nice github.io website__
* xallarap - use case, unit test, and code itself (see below for references)
* binary lens orbital motion - use case, unit test, and code itself
* binary calculations - define if s is relative to total mass etc.
* should BinaryLens() accept source\_x/y as lists or arrays?
* correct JPL Horizons => CSV file format; also example usage
* _t\_* instead of rho_
* _Martin's FSBL code_
* _Cassan 2008 binary lens parameters_
* conversions to BJD from HJD, JD etc. ([astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections))
* for plotting functions option to pass pyplot.Axis and pyplot.Figure instances and call e.g. Axis.scatter() instead of pyplot.scatter()
* subplots with shared X-axis (plt.subplots(2, 1, sharex=True, gridspec\_kw={'height\_ratios': [4, 1]}, figsize=???, dpi=100))
* BJD\_TDB in satellite ephemeris [astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections)
* _faster FSPL with LD_
* function for center of mass shift (currently: shift\_x in trajectory.py, x\_shift in binarylens.py, xcm\_offset in caustics.py)
* single Event can have many instances of Model associated with it
* add unit tests for Horizons and MulensData.satellite\_skycoord
* Caustics.\_calculate - optimize using vectors instead of a loop
* check if Horizons e-mail is for correct satellite
* Are we consistent with PEP8? [check here](http://pep8online.com/)
* use lazy loading in MagnificationCurve.magnification and/or Model.magnification
* _guessing parameters of PSPL model_
* fluxes fixed in chi^2 calculation
* annual parallax calculation - verify with VBBL
* when checking units use Unit.physical\_type
* better import of the module so that all main classes are accessible (use \_\_all\_\_ = [...] in all files?)
* Fit() should use marginalized distributions of fluxes (if those are from linear fits)
* Add \_\_repr\_\_ functions to Lens and Source
* _Model.set\_parameters() should remember previously set values (of course unless they're overwritten)_
* In Lens, add checks for new\_mass as an astropy.units.Quantity and
  use solMass as default if not set.
* FSPL ray shooting (ala getmag\_rs\_single.f)
* t\_eff as a parameter - see [Andy's paper](https://arxiv.org/abs/1312.6692)
* Transform t\_E and other parameters between geocentric and heliocentric frames.
* Errorbar scaling, in particular the two parameter.
* get gamma/u LD coefs from Claret papers etc.
* Class Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
* Class Model should not allow accessing attributes that shouldn't be there, eg., q for single lens case.
* Research decorators (e.g. @property) - Can we have both print(model.s) and print(model.s(time))?
* are we fully ok with astropy license?

### reStructuredText:
[1st tutorial] (http://gisellezeno.com/tutorials/sphinx-for-python-documentation.html)

[2nd tutorial](http://www.sphinx-doc.org/en/stable/rest.html)

[example](https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html)

### Xallarap references:

ob07368 - [Sumi et al. 2010](http://adsabs.harvard.edu/abs/2010ApJ...710.1641S) and [Suzuki et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...833..145S)

ob150845 = mb15277 - Calen leads

mb10328 - [Furusawa et al. 2013](http://adsabs.harvard.edu/abs/2013ApJ...779...91F)

[Poindexter et al. 2005](http://adsabs.harvard.edu/abs/2005ApJ...633..914P) - 23% of events are affected by xallarap

ob9919 - [Smith et al. 2002](http://adsabs.harvard.edu/abs/2002MNRAS.336..670S)

[Han & Gould 1997](http://adsabs.harvard.edu/abs/1997ApJ...480..196H)

[Griest & Hu 1992](http://adsabs.harvard.edu/abs/1992ApJ...397..362G)


