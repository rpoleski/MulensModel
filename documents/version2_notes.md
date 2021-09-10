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

A secondary goal of this version update was to clean up other aspects of the 
code. For example, clarifying keyword and function 
names. We also explicitly created a utils.PlotUtils() so that convenience 
functions for plotting that were used by multiple classes were collected in one
place. Finally, the Event() class no longer stores the best_chi2() or
best_chi2_parameters() from previous updates. These properties are intrisic to
model minimization and not to an Event() object, so should not be properties
of an Event().

# Main Differences

## mm.Fit() vs. mm.FitData()

### mm.Fit()

The old mm.Fit() class acted on all datasets at once and took magnification 
rather than a model as its input. Hence, setting conditions (such as zero 
blending) for individual datasets was not possible with the existing class. In 
addition, because Fit() took magnification as an input, it became a property of
mm.Model() and thus, required that mm.Model() also took datasets and a property.
This did not make logical sense and resulted in circular code.

This class has been DEPRECATED in favor of mm.FitData(). However, if you only
combined datasets with a mm.Model() using the mm.Event() class, you should not
notice much, if any, difference (except increased functionality).

### mm.FitData()

The mm.FitData() class combines a model with a single dataset. This resolves 
the circular logic and allows different fitting conditions (such as zero 
blending) to be applied to one (or several), but not all datasets. Fits for
multiple datasets are combined in the mm.Event() class, i.e., mm.Event() 
contains a new property mm.Event.fits, which is a list of FitData() objects, one
for each dataset.

The blend flux, source flux, and/or the source flux ratio (for two source fits),
can be fixed for an individual FitData() object. For a fit with multiple 
datasets, this information is supplied as a dictionary to the Event() object and 
passed to the relevant FitData() object as necessary.

## mm.Event()

### Fitting API and Fluxes

- Fit() --> FitData()
- ADD: fits
- ADD: fix_blend_flux, fix_source_flux
- ADD: fix_source_flux_ratio
- ADD: fit_fluxes()
- ADD: fluxes, source_fluxes, blend_fluxes
- ADD: get_flux_for_dataset(dataset)

### Transfer Actions that Created Circular Dependencies from mm.Model() to mm.Event()
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

- chi2_gradient() --> calculate_chi2_gradient()
- fit_blending (keyword used by get_chi2(), get_chi2_for_dataset(), 
    get_chi2_per_point(), get_chi2_gradient() ): this should be controlled by fix_blend_flux instead.
- mm.Event().get_ref_fluxes(): data_ref keyword will be deprecated, because
    there is now a get_flux_for_dataset() function.

## mm.MagnificationCurve()

### Remove propery .magnification, which duplicates .get_magnification()

- REMOVE: magnification property

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

### Clarify Keyword/Function Names
- plot_magnification(): 
  - flux_ratio_constraint --> source_flux_ratio
  - DEPRECATED: fit_blending keyword
- plot_lc():
  - f_source --> source_flux
  - f_blend --> blend_flux
  - flux_ratio_constraint --> source_flux_ratio
- magnification() --> get_magnification() and:
  - flux_ratio_constraint --> source_flux_ratio
  - ADD: bandpass

### Other Changes
- REMOVE: set_source_flux_ratio(), set_source_flux_ratio_for_band()

## mm.MulensData():

- plot():
  - model keyword is DEPRECATED.
  
## mm.utils.PlotUtils():

New class of utility functions useful for plotting (collected from elsewhere in 
the code):
- get_y_value_y_err()
- find_subtract()
- find_subtract_xlabel()
- get_color_differences()

## Examples That Reflect These Changes

- example_01_models.py
  - f_source, f_blend --> source_flux, blend_flux
- example_02_fitting.py
  - data_ref replaced by explict FitData instance.
- example_06_fit_parallax_EMCEE.py
  - model.set_datasets instances --> Event() instances
  - model.plot_lc() --> Event.plot_model()
- example_09_gradient_fitting.py
  - event.chi2_gradient() --> event.get_chi2_gradient()
  - model.plot_lc(dataref=data) --> model.plot_lc(source_flux=value, 
    blend_flux=value) where the values are extracted using 
    Event().get_flux_for_dataset(data).
- example_10_fitting_and_fluxes.py
  - event.fit.flux_of_sources(dataset), event.fit.blending_flux(dataset) -->
    event.get_flux_for_dataset(dataset)
- example_11_binary_source.py
  - event.model.set_source_flux_ratio(theta_) --> 
    event.fix_source_flux_ratio ={my_dataset: theta_}
  - model.magnification() --> model.get_magnification()
- example_12_fit_satellite_parallax_EMCEE.py
  - added my_event.fit_fluxes()
  - my_event.model.get_ref_fluxes() --> 
    my_event.get_flux_for_dataset(data_ground)
  - f_source, f_blend --> source_flux, blend_flux
- example_17_1L2S_plotting.py
  - model.magnification() --> model.get_magnification()
- example_18_simulate.py
  - model.magnification() --> model.get_magnification()

Also, because best_chi2 and best_chi2_parameters are no longer properties of 
Event(), we changed how we store and access that information in Examples 06, 10,
 11, 12, 13, 15.

# Other Changes

## New Examples:
- example_17_1L2S_plotting.py
- example_16: added additional yaml files.
