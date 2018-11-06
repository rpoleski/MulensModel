# What we want to change when going from v1.X.Y to v2.0.0?

Once the changes are accepted to be made, mark them in the code using warnings.warn("XXX", FutureWarning) and note it here.

### Major changes:

 * Remove datasets from Model - this moves all plotting functions (except plot\_model) to Event. This causes changes in:
   * Model.data\_magnification()
   * Model.get\_residuals()
   * Model.get\_ref\_fluxes() - possibly
   * Model.reset\_plot\_properties() will be removed
 * fit\_blending maybe should be class property, not an option, or maybe even it should be MulensData property! - This change should simplify Event and Fit functions where we see many "if fit\_blending:"

### Minor changes:
 * Delete ModelParamters.pi\_E and leave pi\_E\_N and pi\_E\_E - it's not really used and just complicates the code inside
 * Remove ModelParameters.as\_dict() because it's the same as ModelParameters.parameters
 * ephemerides\_file -> ephemeris\_file
 * Model.plot\_\* - remove \_list parameters etc.
 * Model.get\_residuals should have keyword phot\_fmt, not type to be consistent with other functions

### Yet unsorted/undecided:
 * check all NotImplementedError and maybe remove some functions/options
 * new class for a collection of datasets to make looping over datasets easier; also there will be data\_ref defined
 * VBBL and AC in subdirs of source/MulensModel/ for easier install
 * MulensData can plot itself - _this seems to be already done_
 * stop updating data\_ref
 * the same order of arguments in plotting functions (if possible)

