### What we need to do for binary sources

* finish use cases (Model.set\_source\_flux\_ratio())
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* ModelParameters - Can child properties be changed?
* Event tries to access self.model.\_fit or one of Model functions passes Fit as additional output - make sure which one exactly
* user can define t\_eff\_1 instead of u\_0\_1 etc.
* Model.magnification() - allow float flux\_ratio\_constraint
* implement Model.same\_source\_flux\_ratio(band="I")
* add binary source parameters to \_valid\_parameters in modelparameters.py and check that which\_parameters() works properly
* check that fit\_blending option is properly considered
* make sure that you can plot model without data (with q\_f specified)
* t\_0\_par and t\_0\_kep - what are default values for binary sources?
* check binary source parallax models
* when finding q\_f via regression, make sure the source fluxes are not negative, at least give warning
* when finging q\_f via regression, how we can get q\_f value?
* ModelParameters.\_\_repr\_\_ has to be updated
* example with q\_f as a MCMC chain parameter
* in modelparameters.py in first 3 functions, what should be names of models (like PSPL, FSBL) for 2 sources where 1 or 2 of them are finite etc.?
* We would like to plot e.g. A\_1/(1+q\_f) or A\_1 etc. - how to access these quantities?
* PEP8 for files with most changes
* go to master branch and review [documents/binary_source_notes.md](https://github.com/rpoleski/MulensModel/blob/master/documents/binary_source_notes.md)

If merging and some parts are not finished, then add ```raise NotImplementedError()``` where needed.

When merging, also update developers_board.md, because some tasks should have been removed from there. Or maybe add something.

