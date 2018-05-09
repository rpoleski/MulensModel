## May goals:
1. binary source use cases
2. errorbar scaling in MulensData
3. finish improved passing of color/label/etc from MulensData to Model
4. website with notes on data challenge and links related to it
5. easier access to all the fluxes from Event
6. improved usage of polynomial root solvers


## Specific tasks to be performed
**boldfaced** tasks are most important because requested by the users

* Install
  * makefile for Windows (basic instructions exist already)
* Documentation
  * Sagan workshop hands-on activity in MM
  * Add \_\_repr\_\_ functions to Lens and Source
  * files not yet well documented: 
    * RA & Dec in coordinates.py (maybe also code it better)
    * coordinates.py
    * utils.py 
    * modelparameters.py
* Effects
  * Finite Source
    * FSPL with low magnification - do [Witt & Mao 94](http://adsabs.harvard.edu/abs/1994ApJ...430..505W) or [Witt 95](http://adsabs.harvard.edu/abs/1995ApJ...449...42W) give the right formulas?
    * faster FSPL with LD
    * FSPL ray shooting (ala getmag\_rs\_single.f)
    * Yoo+04 full formalism 
    * [2] get gamma/u LD coeffs from Claret papers etc.
    * Full formalism of [Lee+09](http://adsabs.harvard.edu/abs/2009ApJ...695..200L)
  * Higher Order Effects
    * xallarap (see below for references)
* Parameterization
  * Cassan 2008 binary lens parameters
  * Albrow et al. 1999 (also Cassan 2008 Sec. 5)
  * t\_eff as a parameter - see [Andy's paper](https://arxiv.org/abs/1312.6692) and maybe also other from [Jen's 2012 paper](http://adsabs.harvard.edu/abs/2012ApJ...755..102Y), i.e., f\_lim=f\_s/u\_0 and q\*t\_E
* Function Improvements/Expansion
  * Binary Lens:
    * should BinaryLens() accept source\_x/y as lists or arrays?
    * function for center of mass shift (currently: shift\_x in trajectory.py, x\_shift in binarylens.py, xcm\_offset in caustics.py)
    * topology of caustics based on (s,q)
    * central and planetary caustic properties: [Chung et al. 2005](http://adsabs.harvard.edu/abs/2005ApJ...630..535C) and [Han 2006](http://adsabs.harvard.edu/abs/2006ApJ...638.1080H)
    * consider using Utils.complex\_fsum() in BinaryLens functions: \_polynomial\_roots\_ok\_WM95() and \_jacobian\_determinant\_ok\_WM95()
  * Caustics.\_calculate - optimize using vectors instead of a loop
  * Caustic calculations using [Erdl & Schneider 1993](http://adsabs.harvard.edu/abs/1993A%26A...268..453E) approach
  * Event class:
    * Event should sync information on which of the 3 types of parallax are used, so that if it's specified for event, then there will be exception if one dataset is missing earth\_coords etc. In general there should be some way to make sure which parallax types are used in which calculation of magnification. 
    * Class Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
    * **Allow fluxes to be fixed in chi^2 calculation (e.g. given a particular fs, fb, which you might want to do if you want fs as a chain parameter); also think how it will work for binary sources**
    * **give access to all fluxes without changing data\_ref**
    * reduce calls to Fit.fit\_fluxes()
    * add finite source in chi2\_gradient()
    * chi2\_gradient() should cope NaN values in a way similar to get\_chi2()
    * check all functions that should pass fit\_blending parameter
    * chi2 with maximum value provided - if the chi2 for point-source gives chi2 larger than specified limit, then finite source calculations are not undertaken (this should significantly speed-up MultiNest)
  * Fit:
    * should use marginalized distributions of fluxes (if those are from linear fits); JCY - it needs UC
  * Horizons:
    * JPL Horizons
      * correct JPL Horizons => CSV file format; also example usage
      * check if Horizons e-mail is for correct satellite
    * BJD
      * conversions to BJD from HJD, JD etc. ([astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections))
      * BJD\_TDB in satellite ephemeris [astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections)
  * Lens:
    * \_\_repr\_\_ function needs work                                         
    * a\_proj, couples with source distance in mulensmodel to determine s.  
    * 2-body example 3 is missing s. Why? Does that work?                  
    * problem with tracking number of masses, esp when successively defining masses (see test\_Lens.py)
    * implement triple+ systems  
  * MagnificationCurve:
    * re-write magnification() to use lazy loading (here or in model.py)
  * Model:  
    * Model.set\_parameters() should remember previously set values (of course unless they're overwritten)
    * Class Model should not allow accessing attributes that shouldn't be there, eg., q for single lens case.
    * Function that prints RA, Dec, t\_0\_par, t\_0\_kep, types of parallaxes turned on, and textual description of type of model
    * plot\_trajectory() should use actual trajectory, not alpha because it may be confusing when orbital motion and parallax are used
    * plot\_trajectory() - mark epochs using colorscale? Maybe it should be passed by kwargs (if so, then add example)
    * Should get\_satellite\_coords() use caching?
    * we should have versions of all plot functions to use magnifications instead of magnitudes; also add access via Event
  * ModelParameters:
    * check that non-existing parameters are not specified e.g. t0
    * check that minimal parameters needed to specify a model are defined
    * Transform t\_E and other parameters between geocentric and heliocentric frames.
    * option to return t\_E, alpha, dalpha\_dt etc. as floats instead of astropy.quantities
  * MulensData:
    * **add label/color/... which is passed to all the matplotlib functions and hence allows to show legend in easy way**
    * **Errorbar scaling, in particular the two parameter.**
    * add version of n\_epochs that uses only good epochs
    * read settings from file header: flux vs. mag, filter, satellite info
  * PointLens:
    * get\_pspl\_magnification() - change it to operate on u^2, not u, so that np.sqrt() calls are reduced
    * 1+2/u^4 approximation for very large u
  * SatelliteSkyCoord:
    * attach magnification\_methods to SatelliteSkyCoord so that they overwrite Model and MagnificationCurve settings when given SatelliteSkyCoord is used
  * Trajectory:
    * **\_get\_delta\_satellite() should be using self.times**
  * Utils:
    * in np.any() ifs give more information in warning e.g., "out of 1234 values provided, the fails are: 12, 345, 678 (0-based)"
    * add u(a) function: u = np.sqrt(2A/np.sqrt(A^2-1.) - 2.)
  * Plotting
    * for plotting functions option to pass pyplot.Axis and pyplot.Figure instances and call e.g. Axis.scatter() instead of pyplot.scatter(); for a simple example see [here](https://github.com/rpoleski/K2-CPM/blob/master/source/K2CPM/plot_utils.py)
    * subplots with shared X-axis (plt.subplots(2, 1, sharex=True, gridspec\_kw={'height\_ratios': [4, 1]}, figsize=???, dpi=100)) - start in Example 5
  * Examples:
    * **chi2 per dataset**
    * **scipy.curve\_fit() and print parameter uncertainties**
    * PSPL fitting with gradient
    * **fs,fb uncertainties**
    * add example that shows 'log\_' in the name of the parameter; central caustic anomaly planet would be best,
    * add illustration on how to remove airmass trends
    * add example of fitting PSPL model using [Albrow (2004)](http://adsabs.harvard.edu/abs/2004ApJ...607..821A) method
    * **corner plots; they require [corner](https://github.com/dfm/corner.py), [pyhdust](https://github.com/danmoser/pyhdust), or [pyGTC](https://pypi.org/project/pyGTC/)**
  * Miscellaneous:
    * when checking units use Unit.physical\_type - search for physical\_type in mulensobjects/lens.py as an example; to find places to be changed search for "isinstance" (to find these places run grep isinstance \*py mulensobjects/\*py | grep Quantity
    * use lazy loading in MagnificationCurve.magnification and/or Model.magnification
    * guessing parameters of PSPL model ([Kim+17](https://arxiv.org/abs/1703.06883) as an example)
    * add calculation of Caustic Region of Influence (CROIN) - [Penny 2014](http://adsabs.harvard.edu/abs/2014ApJ...790..142Y)
    * anything from use cases that does not work yet -- see TODO.md file
    * **plotting data in MulensData (also update PSPL tutorial)**
    * interaction with fitting routines - see [list of them](https://arxiv.org/abs/1711.03329)
    * caching of results in trajectory.py should stop at some point - if the user changes t\_0\_par or coords, then there is no point in remembering huge indexes (whole self.times)
    * profile the code (python -m cProfile script.py)
    * Leap seconds library - [barycorrpy](https://arxiv.org/abs/1801.01634)
* Other Tests:
  * add unit tests for Horizons and MulensData.satellite\_skycoord
  * Coordinates - write tests, possibly remove test\_Coords.py
  * annual parallax calculation - verify with VBBL
  * t\_eff is not tested
* Style/Architecture:
  * Are we consistent with PEP8? [check here](http://pep8online.com/) - last time checked on 28 Feb 2018 (but didn't include tests)
  * better import of the module so that all main classes are accessible (use \_\_all\_\_ = [...] in all files?)
  * Utils - Make subpackage/submodules that group related functions (e.g. flux2mag conversions)?

### reStructuredText:
[1st tutorial] (http://gisellezeno.com/tutorials/sphinx-for-python-documentation.html)

[2nd tutorial](http://www.sphinx-doc.org/en/stable/rest.html)

[example](https://thomas-cokelaer.info/tutorials/sphinx/docstring_python.html)

### Xallarap references:

[Griest & Hu 1992](http://adsabs.harvard.edu/abs/1992ApJ...397..362G),
[Han & Gould 1997](http://adsabs.harvard.edu/abs/1997ApJ...480..196H),
[Dominik 1998](http://adsabs.harvard.edu/abs/1998A%26A...329..361D),
[Ghosh et al. 2004](http://adsabs.harvard.edu/abs/2004ApJ...615..450G),
[Jiang et al. 2004](http://adsabs.harvard.edu/abs/2004ApJ...617.1307J)

ob9919 - [Smith et al. 2002](http://adsabs.harvard.edu/abs/2002MNRAS.336..670S)

[Poindexter et al. 2005](http://adsabs.harvard.edu/abs/2005ApJ...633..914P) - 23% of events are affected by xallarap

ob07368 - [Sumi et al. 2010](http://adsabs.harvard.edu/abs/2010ApJ...710.1641S) and [Suzuki et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...833..145S)

mb10328 - [Furusawa et al. 2013](http://adsabs.harvard.edu/abs/2013ApJ...779...91F)

ob150845 = mb15277 - Calen leads

