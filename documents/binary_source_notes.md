# Action items:


* user can define t\_eff\_1 instead of u\_0\_1 etc. - **easy**
* make sure that you can plot model without data (with q\_f specified) **easy**
* make sure that you can plot model without data and specify q\_f using bands **easy**
* when finding q\_f via regression, make sure the source fluxes are not negative, at least give warning that includes print(model)
* Model.get\_ref\_fluxes - see notes there and decide how many values should be returned
* add fit\_blending parameter for binary sources only: Model.magnification(), data\_magnification(), get\_data\_magnification(), Model.get\_ref\_fluxes(), Model.get\_residuals() AND plot functions; Event.get\_chi2\_for\_dataset() and get\_chi2\_per\_point(), get\_chi2() as well
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly **maybe solved already**
* Event.best\_chi2\_parameters does not remember value set using Event.model.set\_source\_flux\_ratio()
* finish use cases: 21, 26, and 26b (Model.set\_source\_flux\_ratio()) **important**
* add binary source parameters to \_valid\_parameters in modelparameters.py and check that which\_parameters() works properly
* in modelparameters.py in first 3 functions, what should be names of models (like PSPL, FSBL) for 2 sources where 1 or 2 of them are finite etc.?
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* t\_0\_par and t\_0\_kep - Jen suggests t\_0\_1 as most obvious choice
* check 2L2S models [reference](https://ui.adsabs.harvard.edu/abs/2018AJ....155..141B/abstract)
* ModelParameters - How to forbid changing "child" properties in ModelParameters by the user? - re-define dict class?
* when finding q\_f via regression, how we can get q\_f value? **important**
* We would like to plot e.g. A\_1/(1+q\_f) or A\_1 etc. - how to access these quantities?
* Event.chi2\_gradient

**"raise NotImplementedError()" where needed**

# General comments:

* q\_f is needed in Model because without it, there is no effective magnification! Maybe q\_f can be ModelParameters (not Model) property. Looking at current API it seems q\_f is more similar to ModelParameters properties (just floats in most cases) than Model (complicated functions in most cases). On the other hand, q\_f is complicated thing - depends on a dataset or band and requires significant calculations, hence maybe should be in Model.


High level functions that need changes - ALREADY DONE:

* Model:
  * magnification()
  * data\_magnification()
  * get\_data\_magnification()
  * get\_ref\_fluxes()
  * get\_residuals()
  * functions for setting methods (this is less important for now): set\_default\_magnification\_method(), set\_magnification\_methods(), and set\_magnification\_methods\_parameters()
  * all plotting functions - already done: plot\_magnification(), plot\_lc(), plot\_residuals()
* Event:
  * get\_chi2()
  * get\_chi2\_for\_dataset() - **THIS ONE FAILS NOW** (at least with set\_source\_flux\_ratio\_for\_band())
  * get\_chi2\_per\_point()
  * get\_ref\_fluxes()
  * all plotting functions (they call Model functions)

Also note that there is Fit.get\_n\_sources() function that should be taken care off.

Things related to binary source that we'll do in future:

* are satellite data causing any additional problems
* xallarap
* binary-lens/binary-source models
* there may be different limb darkening coeffs for each source - this would affect MulensData
* three sources? We already had to consider them to fully analyze events that turn out to be 2L2S or 3L1S
* different t\_E (and alpha in the case of binary lens) for each source - see [Han+17 on ob160263](https://ui.adsabs.harvard.edu/abs/2017AJ....154..133H/abstract) and [Bennett+18 on mb10117](https://ui.adsabs.harvard.edu/abs/2018AJ....155..141B/abstract)
* it would be nice to have a plotting function that plots combined model, but also A1\*F\_S1+F\_B and A2\*F\_S2+F\_B, so that we can well see the contribution of each source
* if there are different datasets in the same filter, then there is a function that makes flux ratio very similar (or exactly the same) and fitted via regression

