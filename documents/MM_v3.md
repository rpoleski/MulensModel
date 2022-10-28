# What we want to change when going from v2.X.Y to v3.0.0?

Once the changes are accepted to be made, **mark them in the code using warnings.warn("XXX", FutureWarning)** and note it here. Also release a version that differs from previous one only in these warnings - this will allow users to correct their codes.  Also give **suggested changes in warnings**.

### Major changes:

 * search for all "deprecated" are remove it
 * rename Caustics -> CausticsBinary and CausticsWithShear -> CausticsBinaryWithShear (and files) so that they're consistent with CausticsPointWithShear

???

### Minor changes:
 * Delete ModelParameters.pi\_E and leave pi\_E\_N and pi\_E\_E - it is not really used and just complicates the code inside
 * Remove ModelParameters.as\_dict() because it is the same as ModelParameters.parameters
 * `ModelParameters.is_static` -> `is_lens_static`
 * ephemerides\_file -> ephemeris\_file - maybe
 * Model.get\_residuals should have keyword phot\_fmt, not type to be consistent with other functions
 * test\_MulensData.py - in test\_copy() remove warnings.catch\_warnings() because you remove coords, ra, and dec from init

### Yet unsorted/undecided:
 * `Model.set\_times()` - `n\_epochs` should be None as default, so that we can check if both dt and `n\_epochs` were set
 * Caustics.get\_caustics() should return np.arrays, not lists
 * check all NotImplementedError and maybe remove some functions/options
 * somehow change which\_parameters() in modelparameters.py - maybe remove
 * new class for a collection of datasets to make looping over datasets easier; also there will be data\_ref defined
 * the same order of arguments in plotting functions (if possible)
 * ModelParameters - all parameters should be float, not astropy.Quantity objects
 * see (this comment by Jen)[https://github.com/rpoleski/MulensModel/pull/15#issuecomment-1080879537] on how magnification methods are named and called in different parts of the code

### Version 4:
 * Add an Observatory class.
