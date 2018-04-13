* logically, the 2 sources will be specified in ModelParameters
* current UC: have t\_0\_1 and t\_0\_2 instead of t\_0, same for u\_0
* think how we calculate effective magnification - is q\_f ModelParameters property? Probably No, because it can be fitted via regression; maybe the best is to have q\_f as a Model property that can be either fix or fitted

Hence, the best approach I see is to have in Model and ModelParameters private properties that will be instances of these classes but for single sources. Then the main part of Model.magnification will look something like:

```python
mag_0 = magnification_curve_0.magnification
mag_1 = magnification_curve_1.magnification
magnification = (mag_0 * f_s_0 + mag_1 * f_s_1) / (f_s_0 + f_s_1)
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
  * all plotting functions (they call Model functions

Things related to binary source that we'll do in future:

* different methods for both sources - we need to consider binary source with finite source for events that look like lowest mass planets
* is satellite data and binary source causing any additional problems
* xallarap
* binary-lens/binary-source models
* there may be different limb darkening coeffs for each source

