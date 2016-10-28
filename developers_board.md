## High level decisions we should make:
(starting from the most important ones)

1. What conventions do we want for time vector in input files? There are 2 problems currently: 1) astropy supports only JD, not HJD, and 2) WFIRST will be observing post JD=2460000, so most frequently used shorthand format JD' = JD-2450000. = ABCD.XXXXX will have to be modified or extended to 1ABCD.XXXXX. 
1. How to handle changes in origin of the coordinate system? Internally we're working in the center of mass, but fitting is sometimes much faster if the system is relative to, e.g., planetary caustic or a cusp. Also think how this should be done for triple lenses. 
1. How to handle full reparametrization of the model? Cassan 2008 is particular example. 

## October 2016 goals:
1.    Allow fits with fixed no blending.
2.    Make sure get_chi2() works properly for simple binary source model.
3.    Finish MulensData() definition that includes JD-HJD conversion.
4.    Expand unit tests (and code they test) for Event and Model, write them for Fit class.
5.    Calculate source trajectory with annual parallax effect.

## Specific tasks to be performed

* change "type(val) is SomeType" to "isinstance(val, SomeType)"
* one class per file?
* no unit tests for private functions: \_fun()
* long lines & PEP8
* get rid off get_jd_zeropoint from MulensData and its tests
* write tests for MulensData(data_list=...)
* in unit tests if you want to assert that exception was raised then use [these](http://stackoverflow.com/questions/129507/how-do-you-test-that-a-python-function-throws-an-exception) methods

## Use cases to be written 

* Anything from "high level decisions" above.
* Fitting PSPL with free blending and fixed no blending.
* Scaling of observed data to a scale of other dataset. We normally do it to transform follow-up data to survey magnitude scale so that they can be presented on a single plot. 
* Class Model should not allow accesing attributes that shouldn't be there, eg., q for single lens case.
* Transform t_E and other parameters between geocentric and heliocentric frames.
* EMCEE example - see [website](http://dan.iel.fm/emcee/current/user/line/).
* Errorbar scaling, in particular the two parameter.
* Source limb darkening profile: use of gamma and u conventions, obtaining the value from outside sources (Claret papers). 
* RA & Dec in Model - get galactic (.galactic) and elliptical coordinates (from astropy.coordinates import GeocentricTrueEcliptic; .transform_to(GeocentricTrueEcliptic)).


