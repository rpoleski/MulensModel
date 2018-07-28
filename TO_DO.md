### What we need to do for binary sources

* finish use cases
* should we use q\_f or flux\_ratio etc. (what about triple source events?)
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* ModelParameters.n\_sources - make it more intelligent, best if initialized in \_\_init\_\_
* t\_0\_X and u\_0\_X in ModelParameters.\_check\_valid\_combination() - Maybe it should be better to initialize each sub-instance and there will be checks if everything is fine? Otherwise we will have lots of basically copied statements (rho and tstar etc.) - see next task; ModelParameters.\_check\_valid\_combination() should be only checked for single sources? - rename the function
* Model.\_set\_each\_source\_parameters() should probably be moved to modelparameters; Then how to update them whenever any \_1 or \_2 parameter is updated?
* Model.magnification() - better name for source\_fluxes parameter
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly
* Event.get\_chi2\_for\_dataset() - unit test
* rho or t\_star for one or 2 sources
* Model.\_magnification\_2\_sources() - check for fixed q\_f; also implement single q\_f for all datasets provided by the user
* in model.py internal variables: q\_f -> source\_flux\_ratio
* user can define t\_eff\_1 instead of u\_0\_1 etc.
* implement Model.same\_source\_flux\_ratio(band="I")
* add binary source parameters to \_valid\_parameters in modelparameters.py and check that which\_parameters() works properly
* check that fit\_blending option is properly considered
* t\_0\_par and t\_0\_kep - what are default values for binary sources?
* ModelParameters.\_\_repr\_\_ has to be updated
* example with q\_f as a MCMC chain parameter
* in modelparameters.py in first 3 functions, what should be names of models (like PSPL, FSBL) for 2 sources where 1 or 2 of them are finite etc.?
* got to master branch and review [documents/binary_source_notes.md](https://github.com/rpoleski/MulensModel/blob/master/documents/binary_source_notes.md)

If merging and some parts are not finished, then add ```raise NotImplementedError()``` where needed.

When merging, also update developers_board.md, because some tasks should have been removed from there. Or maybe add something.

