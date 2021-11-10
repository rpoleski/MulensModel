# What we want to change when going from v2.X.Y to v3.0.0?

Once the changes are accepted to be made, **mark them in the code using warnings.warn("XXX", FutureWarning)** and note it here. Also release a version that differs from previous one only in these warnings - this will allow users to correct their codes.  Also give **suggested changes in warnings**.

### Major changes:

 * search for all "deprecated" are remove it

???

### Minor changes:
 * Delete ModelParameters.pi\_E and leave pi\_E\_N and pi\_E\_E - it is not really used and just complicates the code inside
 * Remove ModelParameters.as\_dict() because it is the same as ModelParameters.parameters
 * ModelParameters.is_static -> is_lens_static
 * ephemerides\_file -> ephemeris\_file - maybe
 * Model.get\_residuals should have keyword phot\_fmt, not type to be consistent with other functions

### Yet unsorted/undecided:
 * Model.set_times() - n_epochs should be None as default, so that we can check if both dt and n_epochs were set
 * Caustics.get\_caustics() should return np.arrays, not lists
 * check all NotImplementedError and maybe remove some functions/options
 * somehow change which\_parameters() in modelparameters.py - maybe remove
 * new class for a collection of datasets to make looping over datasets easier; also there will be data\_ref defined
 * the same order of arguments in plotting functions (if possible)
 * ModelParameters - all parameters should be float, not astropy.Quantity objects

### Version 4:
 * Add an Observatory class.
