General comments:

* logically, the 2 sources will be specified in ModelParameters
* current UC: have t\_0\_1 and t\_0\_2 instead of t\_0, same for u\_0
* q\_f is needed in Model because without it, there is no effective magnification! Maybe q\_f can be ModelParameters (not Model) property. Looking at current API it seems q\_f is more similar to ModelParameters properties (just floats in most cases) than Model (complicated functions in most cases). 
* q\_f can be fixed or fitted via regression
* Model should somehow remember Fit instance used to find q\_f and then if Event wants to access Fit, then it first checks, if Model has it already. It's better to remember Fit, not just q\_f because it saves least squares function calls.
* We should find better name than q\_f because it's similar to q\_1 that will be used in triple lenses. Maybe source\_flux\_ratio or flux\_ratio?


Hence, the best approach I see is to have in Model and ModelParameters (maybe Event as well) private properties that will be instances of these classes but for single sources. Then the main part of Model.magnification will look something like:

```python
magnification_curve_0 = MagnificationCurve(time, parameters=self.parameters_source_0, ...)
magnification_curve_1 = MagnificationCurve(time, parameters=self.parameters_source_1, ...)
mag_0 = self._magnification_curve_0.magnification
mag_1 = self._magnification_curve_1.magnification
q_f = self._get_q_f(magnification_curve_0, magnification_curve_1) # If fixed, then just 
# returns the value. Otherwise calls function from Event or Fit that finds all 3 fluxes 
# (f_s0, f_s0, and f_b) and returns q_f = f_s0/f_s1.
# If calculates these values, then remember Fit instance, so that Event can access it; think about when to reset this Fit instance.
magnification = (mag_0 + mag_1 * q_f) / (1. + q_f)
return magnification
```

High level functions that need changes:

* Model:
  * magnification()
  * data\_magnification()
  * get\_data\_magnification()
  * get\_ref\_fluxes()
  * get\_residuals()
  * functions for setting methods (this is less important for now): set\_default\_magnification\_method(), set\_magnification\_methods(), and set\_magnification\_methods\_parameters()
  * all plotting functions
* Event:
  * get\_chi2()
  * get\_chi2\_for\_dataset()
  * get\_chi2\_per\_point()
  * get\_ref\_fluxes()
  * all plotting functions (they call Model functions)

Also note that there is Fit.get\_n\_sources() function that should be taken care off.

Things related to binary source that we'll do in future:

* what happens if there are multiple datasets in the same band? In that case, we want to fit the primary source flux and blend flux independently for each dataset, but constrain q\_f to be the same for all datasets with the same band.
* different methods for both sources - we need to consider binary source with finite source for events that look like lowest mass planets; short-term solution can be to change FSPL calculations to PSPL when FS corrections are smaller than 10^-6 or so.
* is satellite data and binary source causing any additional problems
* xallarap
* binary-lens/binary-source models
* there may be different limb darkening coeffs for each source - this would affect MulensData
* chi2 gradient

