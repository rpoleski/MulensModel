### To be discussed:

- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???
- should version be printed on output?
- settings for making plot without a fit - how input should look like for best_model?
- how the output should look like?

### Documentation:
- model['coords']
- fixed parameters
- plotting to screen
- introduce YAML files
- Delta t_0
- binary source

## List of task to be done:

( **boldface** - do this before sending e-mail around)

**NOW - FINISH Mroz+17** - example

- **limit epochs in "best model" plot**
- **methods for Model**
- **some documentation - see above**
- **add one more fitting method?**
- MulensData - just provide *str*
- n_walkers for EMCEE - default is x4 and remove from minimal yaml file
- remove _update_best_model() and extract it from fitting results
- all_parameters in _get_parameters_ordered() and _check_fixed_parameters() - combine in a single one
- note that parameters are re-ordered (maybe in future add option for specifying order)
- datasets - guessing 245/246
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- trace plot
- allow making plots without a fit
- self._plots - check what is there
- Fitting method to be added: scipy.optimize, pymultinest, ???
- allow plotting many random models from posterior
- MulensData() - use try/except with meaningful error message
- plot title
- make plots tighter, i.e., reduce white space
- Add ln_prior values to blob? At some point we will want to save that information in output files
- settings['input_file_root'] = input_file_root - in final function and use it for default output files names
- check if output files (including plots) exists at the begin
- add check if 't_0' is covered by data and give warning if not
- print number of models calculated
- full support of satellite data
- periodic variables - suggest it for alpha, x_caustic_X
- check if data files exist
- allow log10() of parameter
- prior of t_E motivated by Mroz+ papers
- allow turning off flux printing
- warnings on time plotting and data limits - checks for add/subtract 245/246
- example how to run fits on a grid of (s,q)
- allow periodic print of best model etc.
- plot trajectory
- for parallax models check if t_0_par is fixed and give warning, if not
- fits with 0 blending flux for first dataset
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data
- add scipy to _check_imports() - requires siginificant code to be added to _check_imports() in order to find out if t_E prior is used
- if corner could not be imported, then give link to specific file in error message
- flux constraints for binary source models (note that for plotting it is now set to first dataset)
- allow Model.set_magnification_methods_parameters()
