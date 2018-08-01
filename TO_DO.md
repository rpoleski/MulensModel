### What we need to do for binary sources

* finish use cases
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* ModelParameters - Can child properties be changed?
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly
* Model.\_magnification\_2\_sources() - check for fixed q\_f; also implement single q\_f for all datasets provided by the user
* **ModelParameters - define rho\_1 and \_2 using @property**
* ModelParameters.\_check\_valid\_combination() - make sure minimum parameters are defined, also make sure that t\_E is NOT defined both via t\_eff and t\_star
* user can define t\_eff\_1 instead of u\_0\_1 etc.
* implement Model.same\_source\_flux\_ratio(band="I")
* add binary source parameters to \_valid\_parameters in modelparameters.py and check that which\_parameters() works properly
* check that fit\_blending option is properly considered
* t\_0\_par and t\_0\_kep - what are default values for binary sources?
* check binary source parallax models
* ModelParameters.\_\_repr\_\_ has to be updated
* example with q\_f as a MCMC chain parameter
* in modelparameters.py in first 3 functions, what should be names of models (like PSPL, FSBL) for 2 sources where 1 or 2 of them are finite etc.?
* go to master branch and review [documents/binary_source_notes.md](https://github.com/rpoleski/MulensModel/blob/master/documents/binary_source_notes.md)

If merging and some parts are not finished, then add ```raise NotImplementedError()``` where needed.

When merging, also update developers_board.md, because some tasks should have been removed from there. Or maybe add something.

