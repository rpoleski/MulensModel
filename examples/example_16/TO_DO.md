### To be discussed:

- satellite data
- periodic variables
- log10 variables etc.
- flux constraints
- for plots: t_0, \Delta t_0, or t_0 - 2456780 ???
- should version be printed on output?
- settings for making plot without a fit - how input should look like for best_model?
- how the output should look like?

### Documentation:
- Delta t_0
- binary source
- x_caustic_in etc.
- ob08092-o4_prior.yaml

## List of task to be done:

( **boldface** - do this before sending e-mail around)

- **some documentation - see above**
- **add one more fitting method?**
- **corner.py**
- n_walkers for EMCEE - default is x4 and remove from minimal yaml file - "SIMPLIFIES INPUT"
- Mroz+20 - finish
- print fixed parameters at begin or "no fixed parameters", so that full model can be extracted without the input file
- LD coeffs as parameters
- all_parameters in _get_parameters_ordered() and _check_fixed_parameters() - combine in a single one
- note that parameters are re-ordered (maybe in future add option for specifying order)
- datasets - guessing 245/246; plotting as well
- no_negative_blending_flux - only first dataset, or all datasets? Maybe add one more option
- trace plot
- allow plotting multiple models
- for plot script add printing chi2 and fluxes
- allow making plots without a fit
- starting parameters are read from file
- some of the starting values are calculated based on equation given in yaml file, eg. "s: equation 100 / t_E" and then substitute each value of t_E and then use: "exec('out = 100 / 20.12345')" and use variable 'out'; This requires import from math of log, log10, arcsin etc.; make sure "s" in x_caustic_in is not replaced etc.; 
- if Cassan08 paramaterization is used then make sure times are >2450000.
- self._plots - check what is there
- add automatic "obvious" checks on parameters: t_E>0, rho>0, s>0, 1>q>0 - even if they are not provided, then model should be rejected and warning given
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
- Event.get_chi2() - add fit_blending=False option
- allow turning off flux printing
- warnings on time plotting and data limits - checks for add/subtract 245/246
- if code fails during fitting, then it should still print the best model found so far - add try/except in _run_fit()
- example how to run fits on a grid of (s,q)
- allow periodic print of best model etc.
- plot trajectory
- for parallax models check if t_0_par is fixed and give warning, if not
- fits with 0 blending flux for first dataset
- when plotting best model, plot ~100 points based on t_E etc. + all visible epochs in data so that anomalies are not missed etc.
- add scipy to _check_imports() - requires siginificant code to be added to _check_imports() in order to find out if t_E prior is used
- if corner could not be imported, then give link to specific file in error message
- flux constraints for binary source models (note that for plotting it is now set to first dataset)
- method to be used: https://lmfit.github.io/lmfit-py/
- allow Model.set_magnification_methods_parameters()
- add source parameter for Model.set_magnification_methods()
- methods - if only single string is provided, then this is a default method
- print all models
- print current best model - each minute, each nth model etc.
- print every n-th model
