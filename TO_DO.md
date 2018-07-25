### What we need to do for binary sources

* finish use cases
* do we use t\_0\_1 or t\_0\_\_1?
* should we use q\_f or flux\_ratio etc. (what about triple source events?)
* should q\_f calculation depend on MulensData or on band?
* Decide on how to request single source output in double source models and then use it consistently in all plotting functions, magnification functions etc.
* t\_0\_X and u\_0\_X in ModelParameters.\_check\_valid\_combination()
* chi2 and multiple datasets
* Event tries to access self.model.\_fit
* Model.\_set\_each\_source\_parameters() should probably be moved to modelparameters
* Model.\_magnification\_2\_sources() - check for fixed q\_f
* check that fit\_blending option is properly considered
* got to master branch and review [documents/binary_source_notes.md](https://github.com/rpoleski/MulensModel/blob/master/documents/binary_source_notes.md)

