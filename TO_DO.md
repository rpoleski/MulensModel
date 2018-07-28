### What we need to do for binary sources

* finish use cases
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* public access to \_source\_1\_parameters and \_2\_ in ModelParameters and use it in Model.magnification()
* remove Model.\_set\_each\_source\_parameters
* ModelParameters - update self.\_source\_1\_parameters and \_2\_ when parameters are changed for parent; Can child properties be changed?
* Model.magnification() - better name for source\_fluxes parameter
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly
* Event.get\_chi2\_for\_dataset() - unit test
* rho or t\_star for one or 2 sources
* Model.\_magnification\_2\_sources() - check for fixed q\_f; also implement single q\_f for all datasets provided by the user
* in model.py internal variables: q\_f -> source\_flux\_ratio
* ModelParameters.\_check\_valid\_combination() - make sure minimum parameters are defined
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

