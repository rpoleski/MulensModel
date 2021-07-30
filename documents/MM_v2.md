# What we want to change when going from v1.X.Y to v2.0.0?

Once the changes are accepted to be made, **mark them in the code using warnings.warn("XXX", FutureWarning)** and note it here. Also release a version that differs from previous one only in these warnings - this will allow users to correct their codes.  Also give **suggested changes in warnings**.

### Major changes:

 * Remove datasets from Model - this moves all plotting functions (except plot\_model) to Event. Not sure if we really want that - it's currently used for initializing Fit. We have to somehow do lazy loading of magnifications. If so, then his causes changes in:
   * Model.data\_magnification()
   * Model.get\_residuals()
   * Model.get\_ref\_fluxes() - possibly
   * Model.reset\_plot\_properties() will be removed
   * Model.set\_source\_flux\_ratio\_for\_band
 * fit\_blending maybe should be class property, not an option, or maybe even it should be MulensData property! - This change should simplify Event, Model, and Fit functions where we see many "if fit\_blending:" or passing fit\_blending to other functions.
 * example_06: refactor to create a MySampler(emcee.EnsembleSampler) class with functions for extracting median/confidence intervals and best sample, etc.

### Minor changes:
 * Delete ModelParameters.pi\_E and leave pi\_E\_N and pi\_E\_E - it's not really used and just complicates the code inside
 * Remove ModelParameters.as\_dict() because it's the same as ModelParameters.parameters
 * ephemerides\_file -> ephemeris\_file - maybe
 * Model.plot\_\* - remove \_list parameters etc.
 * Model.get\_residuals should have keyword phot\_fmt, not type to be consistent with other functions
 * Model.set\_datasets and .datasets - make sure it's always a list even if only single MulensData object is provided

### Yet unsorted/undecided:
 * Model.set_times() - n_epochs should be None as default, so that we can check if both dt and n_epochs were set
 * Caustics.get\_caustics() should return np.arrays, not lists
 * check all NotImplementedError and maybe remove some functions/options
 * somehow change which\_parameters() in modelparameters.py - maybe remove
 * new class for a collection of datasets to make looping over datasets easier; also there will be data\_ref defined
 * VBBL and AC in subdirs of source/MulensModel/ for easier install ?
 * stop updating data\_ref in Model
 * the same order of arguments in plotting functions (if possible)
 * Model.magnification has parameter gamma - shouldn't it be bandpass?
 * Model.magnification has parameter flux\_ratio\_constraint but there is also Model.set\_source\_flux\_ratio() - probably there should be just one version of that - the latter seems more logical
 * binary source models - currently, I think, the source flux is visible as a 1D array, but its size is (1,) even for binary source events - make sure what's the status and either add full support of (2,) size or just change it to float
 * ModelParameters or Model internally split into classes like: PSPL, BSPL, BSBL and use inheritance; but should those also have parallax information? finite source vs. point source?; We may end up with huge number of classes, which is not good if current solution works well. The only problem is that modelparameters.py is very long. We should also check if binary-lens binary-source models work well.
 * BinaryLens - change all *_magnification() to get_...()
 * ModelParameters - all parameters should be float, not astropy.Quantity objects
 * Event.get_ref_fluxes() and similar functions should not update data_ref ?
 * Event.sum_function - change default to 'numpy.sum'

### Version 3:
 * Add an Observatory class.
