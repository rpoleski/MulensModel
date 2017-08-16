## May 2017 goals:
1. PSPL fitting manual,
2. Limb darkening coefficients -- first we need use case,
4. Resolve parallax issue.


## Specific tasks to be performed
(__boldfaced__ correspond to this month goals)

* __PSPL manual__
* binary calculations - define if s is relative to total mass etc.
* __should BinaryLens() accept source\_x/y as lists or arrays?__
* __correct vector\_tau = (...) in Model.\_trajectory__
* example usage of JPL Horizons
* single Event can have many instances of Model associated with it
* __horizons.py - change "Time(...).jd - 2450000" to something normal__
* __MulensData - change "self.jd - 2450000." to something normal__
* add unit tests for Horizons and MulensData.satellite\_skycoord
* check if Horizons e-mail is for correct satellite
* guessing parameters of PSPL model
* anything from use cases at the end of this page
* fluxes fixed in chi^2 calculation
* annual parallax calculation - verify with VBB
* when checking units use Unit.physical\_type
* one test file per class
* pass datasets from Event to Model or vice versa
* does every function have a description? 
* better import of the module so that all main classes are accesiable
* Fit() should use marginalized distributions of fluxes (if those are from linear fits)
* Reconsider implementation of plotting in use case 08 (perhaps more like use case 02 or based on use case 10)
* make sure Event.\_\_init\_\_ is correct
* Add __repr__ functions to Lens and Source
* In Lens, add checks for new\_mass as an astropy.units.Quantity and
  use solMass as default if not set.
* improve accuracy of test\_annual\_parallax\_calculation() in tests/test\_Model.py
* Classes Model and Event should have not only set\_datasets() methods but also add\_datasets(), i.e. a similar method that appends datasets to self.\_datasets.
* on-line access to JPL Horizons 
* Can model generate magnifications without data?

* Check code with spyder

### Non-functional elements of use cases:
* 01: model does not support ~~time~~, caustics, trajectory
* ~~03: Model does not support source and blend fluxes or magnitudes. Consequently, it also does not support plotting time vs. magnitude.~~
* 03: Model does not support any plotting. Needs ~~plot_lightcurve~~, plot_caustics, plot_trajectory. plot_caustics is actually meant to be a function in Lens.
* 04: definition of ra and dec not supported. see also point above and use case 13.
* 04: Event does not support append (as in append a new dataset)
* 05: Need to download (and implement) default ephemrides for K2
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
* 17: Exposure time, integration over exposure time, time defined as center of integration

## Decisions we should make:

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

