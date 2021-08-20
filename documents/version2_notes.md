# Summary

The major goal of this version update was to allow control of fits for 
individual datasets with the same model. For example, the ability to fix the 
blend flux = 0 for a single dataset but allow it to be a free parameter for
other datasets. Achieving this goal required a significant change in the overall
architecture. In particular, we moved fitting of datasets from inside the 
Model() class to the Event() class. In retrospect, the new architecture makes
more sense because it removes circular dependencies from the Model() class and
makes it completely independent from the data. The basic user who primarily 
accessed fit results through the Event() class should not notice these changes.

A secondary goal of this version update was to clarify keyword and function 
names. 

Finally, the Event() class no longer stores the best_chi2() or
best_chi2_parameters() from previous updates. These properties are intrisic to
model minimization and not to an Event() object, so should not be properties
of an Event().

# Differences

## mm.FitData()

New. Handles each dataset independently. (will write more later.)

## mm.Event()

### Fitting API and Fluxes

- Fit() --> FitData()
- ADD: fits
- ADD: fix_blend_flux, fix_source_flux
- ADD: fix_source_flux_ratio
- ADD: fit_fluxes()
- ADD: fluxes, source_fluxes, blend_fluxes
- ADD: get_flux_for_dataset(dataset)

### Remove Circular Dependencies from mm.Model() 
(and cleanup other plotting functions)

- ADD: utils.PlotUtils()
- ADD: data_ref (optional, defaults to first dataset)
- MOVED contents of many plotting functions from mm.Model() --> mm.Event():
  - plot_model()
  - plot_data()
  - plot_residuals()
  - plot_source_for_datasets()
- MOVED: mm.Model._set_default_colors() --> mm.Event()

### Remove Minimization Properties from mm.Event()

- REMOVE: reset_best_chi2()
- REMOVE: best_chi2, best_chi2_parameters 

### Other Changes

- chi2_gradient() --> calc_chi2_gradient()
- fit_blending (keyword used by get_chi2(), get_chi2_for_dataset(), 
    get_chi2_per_point(), get_chi2_gradient() ): this should be controlled by fix_blend_flux instead.
- mm.Event().get_ref_fluxes(): data_ref keyword will be deprecated, because
    there is now a get_flux_for_dataset() function.
    
## mm.Model()

### Remove Circular Dependencies from mm.Model() 
Combining datasets with a model and all fitting are now handled by mm.Event()

- REMOVE: data_ref keyword.
- REMOVE: fit, source_flux_ratio_constraint, datasets
- plot_trajectory():
  - DEPRECATED: show_data keyword.
- REMOVE: plot_source_for_datasets()
- REMOVE: get_ref_fluxes()
- REMOVE: data_magnification, get_data_magnification(), get_residuals()
- REMOVE: reset_plot_properties(), plot_data(), plot_residuals()
- REMOVE: datasets, set_datasets()

### Clarify Keyword Names
- plot_magnification(): 
  - flux_ratio_constraint --> source_flux_ratio
  - DEPRECATED: fit_blending keyword
- plot_lc():
  - f_source --> source_flux
  - f_blend --> blend_flux
  - flux_ratio_constraint --> source_flux_ratio
- magnification():
  - flux_ratio_constraint --> source_flux_ratio
  - ADD: bandpass

### Other Changes
- REMOVE: set_source_flux_ratio(), set_source_flux_ratio_for_band()