## current goals:
1. ???
2. ???
3. ???

## Specific tasks to be performed
(__boldfaced__ correspond to this month goals; try to put important stuff at the top)

* __save release 0.1.0 for code review - see [github releases](https://help.github.com/articles/creating-releases/)__
* __PSPL manual__
* correct documentation for existing code - start with Model, MulensData, and Event 
* anything from use cases that does not work yet -- see TODO.md file
* copy documentation to some website
* xallarap - use case, unit test, and code itself (see below for references)
* binary lens orbital motion - use case, unit test, and code itself
* binary calculations - define if s is relative to total mass etc.
* should BinaryLens() accept source\_x/y as lists or arrays?
* correct JPL Horizons => CSV file format
* Martin's FSBL code
* convertions to BJD from HJD, JD etc. ([astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections))
* for plotting functions option to pass pyplot.Axis and pyplot.Figure instances and call e.g. Axis.scatter() instead of pyplot.scatter()
* subplots with shared X-axis (plt.subplots(2, 1, sharex=True, gridspec\_kw={'height\_ratios': [4, 1]}, figsize=???, dpi=100))
* BJD\_TDB in satellite ephemeris [astropy link](http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections)
* faster FSPL with LD
* example usage of JPL Horizons
* function for center of mass shift (currently: shift\_x in trajectory.py, x\_shift in binarylens.py, xcm\_offset in caustics.py)
* single Event can have many instances of Model associated with it
* add unit tests for Horizons and MulensData.satellite\_skycoord
* Caustics.\_calculate - optimize using vectors instead of a loop
* check if Horizons e-mail is for correct satellite
* use lazy loading in MagnificationCurve.magnification and/or Model.magnification
* guessing parameters of PSPL model
* fluxes fixed in chi^2 calculation
* annual parallax calculation - verify with VBBL
* when checking units use Unit.physical\_type
* better import of the module so that all main classes are accesible
* Fit() should use marginalized distributions of fluxes (if those are from linear fits)
* Add \_\_repr\_\_ functions to Lens and Source
* In Lens, add checks for new\_mass as an astropy.units.Quantity and
  use solMass as default if not set.
* FSPL ray shooting (ala getmag\_rs\_single.f)
* Class Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
* Research decorators (e.g. @property) - Can we have both print(model.s) and print(model.s(time))?

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

## Decisions we should make:

1. How to handle changes in origin of the coordinate system? Internally we're working in the center of mass, but fitting is sometimes much faster if the system is relative to, e.g., planetary caustic or a cusp. Also think how this should be done for triple lenses. 
1. How to handle full reparametrization of the model? Cassan 2008 is particular example. 


## Use cases to be written 

* Anything from "high level decisions" above.
* Specify that the source has 2 components.
* Scaling of observed data to a scale of other dataset. We normally do it to transform follow-up data to survey magnitude scale so that they can be presented on a single plot. 
* Class Model should not allow accessing attributes that shouldn't be there, eg., q for single lens case.
* Transform t\_E and other parameters between geocentric and heliocentric frames.
* Errorbar scaling, in particular the two parameter.
* Source limb darkening profile: use of gamma and u conventions, obtaining the value from outside sources (Claret papers). 

