## November 2016 goals:
1.    Verify that annual parallax works fine.
2.    Add satellite parallax calculations.
3.    Limit list of specific tasks below to around 10 (or less).

## Specific tasks to be performed
(__boldfaced__ correspond to this month goals)

* __annual parallax calculation - test accuracy__
* __satellite parallax__
* replace Model.reset\_magnification() with remembering parameters for calculated model in an instance of ModelParameters
* when checking units use Unit.physical\_type
* full test of HJD-JD correction in MulensData and Fit
* one test file per class
* pass datasets from Event to Model or vice versa
* check longest files - does every function have a description?
* better import of the module so that all main classes are accesiable
* Fit() should use marginalized distributions of fluxes
* use case 16 - code all coords features
* in Model: Fix ra and dec setters
* Reconsider implementation of plotting in use case 08 (perhaps more
  like use case 02 or based on use_case 10)
* make sure Event.\_\_init\_\_ is correct
* improve accuracy of test\_annual\_parallax\_calculation() in tests/test\_Model.py

### Non-functional elements of use cases:
* 01: Model does not allow parameters to be set in this way (see also use cases 08, 10, 12, 13, 16)
* 02: Entire use case not implemented
  * 02: Model does not support defining the model by the lens and source (see also use case 03).
  * 02: Model does not support plot_lightcurve()
  * 02: Model does not support theta_E (see also use case 03)
* 03: Model does not support source and blend fluxes or magnitudes. Consequently, it also does not support plotting time vs. magnitude.
* 03: Model does not support any plotting. Needs plot_lightcurve, plot_caustics, plot_trajectory. plot_caustics is actually meant to be a function in Lens.
* 04: definition of ra and dec not supported. see also point above and use case 13.
* 04: Event does not support append (as in append a new dataset)
* 05: MulensData does not support satellite
* 06: WFIRST data not implemented, MulensData also does not support bandpass
* 07: This version of defining a Lens is not implemented.
* 08: see above
* 09: Finite source effects not implemented.
* 10: source_flux and blend_flux not supported by Model. I'm not sure
  the plotting would work either.
* 11: Entire use case not implemented
* 12: Event does not support chi2_0
* 13: ModelParameters does not support frame_origin 
* 16: MulensData does not support ra, dec. Event does not support coords (or ra, dec).


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

