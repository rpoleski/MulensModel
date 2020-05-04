### To be discussed:

- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???
- should version be printed on output?
- settings for making plot without a fit - how input should look like for best_model?
- how the output should look like?


## List of task to be done:

( **boldface** - do this before sending e-mail around)

- **coordinates - to allow parallax fits**
- MulensData - just provide *str*
- remove _update_best_model() and extract it from fitting results
- binary source and number of fluxes returned - see _return_ln_prob()
- all_parameters in _get_parameters_ordered() and _check_fixed_parameters() - combine in a single one
- note that parameters are re-ordered (maybe in future add option for specifying order)
- datasets - guessing 245/246
- **data here for tests, so that PATH doesnt change**
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- **limit epochs in "best model" plot**
- **limit Y axis in "best model" plot**
- trace plot
- **methods for Model**
- allow making plots without a fit
- self._plots - check what is there
- Fitting method to be added: scipy.optimize, pymultinest, ???
- allow plotting many random models from posterior
- document plotting to screen
- MulensData() - use try/except with meaningful error message
- plot title
- make plots tighter, i.e., reduce white space
- Add ln_prior values to blob? At some point we will want to save that information in output files
- extend documentation
- settings['input_file_root'] = input_file_root - in final function and use it for default output files names
- check if output files (including plots) exists at the begin
- add check if 't_0' is covered by data and give warning if not
- print number of models calculated
- full support of satellite data
- periodic variables - suggest it for alpha, x_caustic_X
- check if data files exist
- explain Delta t_0 in documentation
- allow log10() of parameter
- prior of t_E motivated by Mroz+ papers
- allow turning off flux printing
- warnings on time plotting and data limits - checks for add/subtract 245/246
- example how to run fits on a grid of (s,q)
- allow periodic print of best model etc.
- fits with 0 blending flux for first dataset
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data
