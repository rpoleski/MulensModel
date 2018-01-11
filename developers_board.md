## Jan goals:
1. prepare v1.0.0 release and A&C paper
2. instantaneous orbital motion of binary lens,
3. Jacobian for PSPL
4. finish new implementation of ModelParameters
5. PointLens class


## Specific tasks to be performed
(__boldfaced__ correspond to this month goals; _italics_ mark task useful for data challenge; try to put important stuff at the top)

[Why?] = JCY thinks this needs a use case

[1] = Necessary for data challenge

[2] = Nice for data challenge

* Install
  * [2] makefile for Windows (basic instructions exist already)
* Documentation
  * Sagan workshop hands-on activity in MM
  * binary calculations
    1. Confirm s is relative to Einstein radius of total mass of the binary lens.
    2. Update documentation to reflect this.
    3. Make code changes as necessary.

  * Add \_\_repr\_\_ functions to Lens and Source
  * [2] files not yet well documented: 
    * RA & Dec in coordinates.py (maybe also code it better)
    * [2] magnificationcurve.py - **needs a list of allowed magnification methods**
    * utils.py 
    * [1] modelparameters.py - **needs more work**
* Effects
  * Finite Source
    * FSPL with low magnification - do [Witt & Mao 94](http://adsabs.harvard.edu/abs/1994ApJ...430..505W) or [Witt 95](http://adsabs.harvard.edu/abs/1995ApJ...449...42W) give the right formulas?
    * [2] _faster FSPL with LD_
    * FSPL ray shooting (ala getmag\_rs\_single.f)
    * [2] get gamma/u LD coefs from Claret papers etc.
  * Higher Order Effects
    * xallarap (see below for references)
      - use case, 
      - unit test, 
      - and code itself 
    * [1] __binary lens orbital motion__
      - use case, 
      - unit test, 
      - and code itself
      * [1] Research decorators (e.g. @property) - Can we have both print(model.s) and print(model.s(time))?
* Parameterization
  * __t\_* instead of rho__
  * _Cassan 2008 binary lens parameters_
  * [1] _dA/dparam for point lens models_
    * ~~[1] use case ~~ --> UC24
  * t\_eff as a parameter - see [Andy's paper](https://arxiv.org/abs/1312.6692) and maybe also other from [Jen's 2012 paper](http://adsabs.harvard.edu/abs/2012ApJ...755..102Y)
* Function Improvements/Expansion
  * Binary Lens:
    * should BinaryLens() accept source\_x/y as lists or arrays?
    * [2] function for center of mass shift (currently: shift\_x in trajectory.py, x\_shift in binarylens.py, xcm\_offset in caustics.py)
  * Caustics.\_calculate - optimize using vectors instead of a loop
  * Event class:
    * Event should sync information on which of the 3 types of parallax are used, so that if it's specified for event, then there will be exception if one dataset is missing earth\_corrds etc. In general there should be some way to make sure which parallax types are used in which calculation of magnification. 
    * Class Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
    * Allow fluxes to be fixed in chi^2 calculation (e.g. given a particular fs, fb, which you might want to do if you want fs as a chain parameter)
  * [Why?] Fit() should use marginalized distributions of fluxes (if those are from linear fits)
  * Horizons:
    * JPL Horizons
      * correct JPL Horizons => CSV file format; also example usage
      * check if Horizons e-mail is for correct satellite
    * BJD
      * conversions to BJD from HJD, JD etc. ([astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections))
      * BJD\_TDB in satellite ephemeris [astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections)
  * Lens:
    * see "To be done:" in mulensobjects/lens.py::Lens docstring
  * Model:  
    * _Model.set\_parameters() should remember previously set values (of course unless they're overwritten)_
    * Class Model should not allow accessing attributes that shouldn't be there, eg., q for single lens case.
  * ModelParameters:
    * should we set t\_E as a Astropy.quantity and then expect to get float?
    * Transform t\_E and other parameters between geocentric and heliocentric frames.
  * Plotting
    * for plotting functions option to pass pyplot.Axis and pyplot.Figure instances and call e.g. Axis.scatter() instead of pyplot.scatter(); for a simple example see [here](https://github.com/rpoleski/K2-CPM/blob/master/source/K2CPM/plot_utils.py)
   * [2] subplots with shared X-axis (plt.subplots(2, 1, sharex=True, gridspec\_kw={'height\_ratios': [4, 1]}, figsize=???, dpi=100))
  * Miscellaneous:
    * when checking units use Unit.physical\_type - search for physical\_type in mulensobjects/lens.py as an example; to find places to be changed search for "isinstance" (to find these places run grep isinstance \*py mulensobjects/\*py | grep Quantity
    * use lazy loading in MagnificationCurve.magnification and/or Model.magnification
    * guessing parameters of PSPL model ([Kim+17](https://arxiv.org/abs/1703.06883) as an example)
    * Errorbar scaling, in particular the two parameter.
    * add calculation of Caustic Region of Influence (CROIN) - [Penny 2014](http://adsabs.harvard.edu/abs/2014ApJ...790..142Y)
    * anything from use cases that does not work yet -- see TODO.md file
    * [2] plotting data in MulensData (also update PSPL tutorial)
    * interaction with fitting routines - see list of them e.g. https://arxiv.org/abs/1711.03329
    * caching of results in trajectory.py should stop at some point - if the user changes t\_0\_par or coords, there there is no point in remembering huge indexes (whole self.times)
    * profile the code (python -m cProfile script.py)
    * Leap seconds library - [barycorrpy](https://arxiv.org/abs/1801.01634)
* Other Tests:
  * add unit tests for Horizons and MulensData.satellite\_skycoord
  * annual parallax calculation - verify with VBBL
* Style/Architecture:
  * Are we consistent with PEP8? [check here](http://pep8online.com/) - last time checked in mid-Dec 2017
  * better import of the module so that all main classes are accessible (use \_\_all\_\_ = [...] in all files?)
* [2] submit to Astronomy&Computing

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


