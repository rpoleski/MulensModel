### What we need to do for binary sources

* finish use cases
* should we use q\_f or flux\_ratio etc. (what about triple source events?)
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* t\_0\_X and u\_0\_X in ModelParameters.\_check\_valid\_combination()
* Model.magnification() - better name for source\_fluxes parameter
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly
* Model.\_set\_each\_source\_parameters() should probably be moved to modelparameters; Then how to update them whenever any \_1 or \_2 parameter is updated?
* Model.\_magnification\_2\_sources() - check for fixed q\_f; also implement single q\_f for all datasets provided by the user
* in model.py internal variables: q\_f -> source\_flux\_ratio
* user can define t\_eff\_1 instead of u\_0\_1 etc.
* implement Model.same\_source\_flux\_ratio(band="I")
* add binary source parameters to \_valid\_parameters in modelparameters.py and check that which\_parameters() works properly
* check that fit\_blending option is properly considered
* ModelParameters.\_\_repr\_\_ has to be updated
* example with q\_f as a MCMC chain parameter
* got to master branch and review [documents/binary_source_notes.md](https://github.com/rpoleski/MulensModel/blob/master/documents/binary_source_notes.md)

If merging and some parts are not finished, then add ```raise NotImplementedError()``` where needed.

