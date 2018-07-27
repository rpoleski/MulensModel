General comments:

* logically, the 2 sources will be specified in ModelParameters (and hence available in Model)
* current UC: have t\_0\_1 and t\_0\_2 instead of t\_0, same for u\_0; 
* q\_f should be Model or ModelParameters property; it seems Model is better, because ModelParameters does not know anything about datasets.
* q\_f can be fixed or fitted via regression (to a single or more datasets)
* Model should remember Fit instance used to find q\_f and then if Event wants to access Fit, then it first checks, if Model has it already. It's better to remember Fit, not just q\_f because it saves least squares function calls.
* Model.magnification() does not currently know which dataset it is using to calculate q\_f via regression - new parameter is needed. 
* More general - q\_f calculation should depend on MulensData or all instances of MulensData in given filter. The problem with MulensData is that if a given dataset does not cover one of the peaks, then flux fit is ill-conditioned (maybe the best solution is to check for it in Fit and possibly rise Exception). There will be option like ```Model.same_source_flux_ratio(band="I")```.
* use source\_flux\_ratio instead of q\_f
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.


Hence, the best approach I see is to have in Model and ModelParameters (maybe Event as well) private properties that will be instances of these classes but for single sources. Then the main part of Model.magnification will look something like:

```python
self._magnification_curve_0 = MagnificationCurve(time, parameters=self.parameters_source_0, ...)
self._magnification_curve_1 = MagnificationCurve(time, parameters=self.parameters_source_1, ...)
mag_0 = self._magnification_curve_0.magnification
mag_1 = self._magnification_curve_1.magnification
if 'q_f' in self._parameters.parameters:
    q_f = self._parameters.q_f
else:
    self._fit = Fit(data=self._datasets, magnification=[mag_0, mag_1]) 
    # Event will try to use this instance of Fit before creating a new one.
    self._fit.fit_fluxes()
    f_s = self._fit.flux_of_sources(SOME_DATASET) # self.datasets[0] as default
    q_f = f_s[1] / f_s[0] 
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
* are satellite data causing any additional problems
* xallarap
* binary-lens/binary-source models
* there may be different limb darkening coeffs for each source - this would affect MulensData
* chi2 gradient
* three sources? We already had to consider them to fully analyze events that turn out to be 2L2S or 3L1S
* different t\_E (and alpha in the case of binary lens) for each source - see [Han+17 on ob160263](http://adsabs.harvard.edu/abs/2017AJ....154..133H) and [Bennett+18 on mb10117](http://adsabs.harvard.edu/abs/2018AJ....155..141B)
* it would be nice to have a plotting function that plots combined model, but also A1\*F\_S1+F\_B and A2\*F\_S2+F\_B, so that we can well see the contribution of each source

