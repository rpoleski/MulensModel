## November 2016 goals:
1.    Verify that annual parallax works fine.
2.    Add satellite parallax calculations.
3.    Limit list of specific tasks below to around 10 (or less).

## Specific tasks to be performed
(__boldfaced__ correspond to this month goals)

* __annual parallax calculation - test accuracy__
* __satellite parallax__
* t\_0\_par is MulensTime instance
* one class per file
* one test file per class
* add check on astropy version minimum 1.2 in MulensData
* write @property for Model that returns Galactic and ecliptic coordinates based on \_coords
* MulensTime._get_date_zeropoint() and MulensData._get_date_zeropoint() - use dictionary 
* remove EMPTY files
* pass datasets from Event to Model or vice versa
* check longest files - does every function have a description?
* add Reduced JD to accepted time formats
* add a check (and warning if found) that data specified are before 1992 or after 2050
* better import of the module so that all main classes are accesiable
* no unit tests for private functions: \_fun()
* Fit() should use marginalized distributions of fluxes
* in unit tests if you want to assert that exception was raised then use [these](http://stackoverflow.com/questions/129507/how-do-you-test-that-a-python-function-throws-an-exception) methods
* all comments are sentences
* use case 16 - code all coords features
* event coords - Galactic and Ecliptic frames
* add Reduced JD and RHJD
* in Model: Fix ra and dec setters
* in Model: allow parameters to be set as in use case 01
* Reconsider implementation of plotting in use case 08

## Decisions we should make:

1. What conventions do we want for time vector in input files? There are 2 problems currently: 1) astropy supports only JD, not HJD, and 2) WFIRST will be observing post JD=2460000, so most frequently used shorthand format JD' = JD-2450000. = ABCD.XXXXX will have to be modified or extended to 1ABCD.XXXXX. 
1. How to handle changes in origin of the coordinate system? Internally we're working in the center of mass, but fitting is sometimes much faster if the system is relative to, e.g., planetary caustic or a cusp. Also think how this should be done for triple lenses. 
1. How to handle full reparametrization of the model? Cassan 2008 is particular example. 

## Use cases to be written 

* Anything from "high level decisions" above.
* Specify that the source has 2 components.
* Scaling of observed data to a scale of other dataset. We normally do it to transform follow-up data to survey magnitude scale so that they can be presented on a single plot. 
* Class Model should not allow accesing attributes that shouldn't be there, eg., q for single lens case.
* Transform t_E and other parameters between geocentric and heliocentric frames.
* Errorbar scaling, in particular the two parameter.
* Source limb darkening profile: use of gamma and u conventions, obtaining the value from outside sources (Claret papers). 
* RA & Dec in Model - get galactic (.galactic) and elliptical coordinates (from astropy.coordinates import GeocentricTrueEcliptic; .transform_to(GeocentricTrueEcliptic)).


